export rising, U, U_dz, U_da, U_dzda

"""
    rising(x, n)

Compute the rising factorial ``(x)ₙ = x(x + 1)⋯(x + n - 1)``.
"""
rising(x::Arb, n::Arb) = Arblib.rising!(zero(x), x, n)
rising(x::Arb, n::Integer) = Arblib.rising!(zero(x), x, convert(UInt, n))
rising(x::ArbSeries, n::Integer) =
    Arblib.rising_ui_series!(zero(x), x, convert(UInt, n), length(x))
rising(x::Acb, n::Acb) = Arblib.rising!(zero(x), x, n)
rising(x::Acb, n::Integer) = Arblib.rising!(zero(x), x, convert(UInt, n))
rising(x::AcbSeries, n::Integer) =
    Arblib.rising_ui_series!(zero(x), x, convert(UInt, n), length(x))
rising(x::T, n::T) where {T<:Union{Float64,ComplexF64}} = Arblib.fpwrap_rising(x, n)
rising(x, n::Integer) = prod(i -> x + i, 0:n-1, init = one(x))

"""
    U(a, b, z)

Compute ``U(a, b, z)``, the confluent hypergeometric function of
the second kind.
"""
U(a, b, z) = U(promote(a, b, z)...)

U(a::Arblib.ArbOrRef, b::Arblib.ArbOrRef, z::Arblib.ArbOrRef) =
    Arblib.hypgeom_u!(zero(z), a, b, z)

U(a::Arblib.AcbOrRef, b::Arblib.AcbOrRef, z::Arblib.AcbOrRef) =
    Arblib.hypgeom_u!(zero(z), a, b, z)

function U(a, b, z::Union{ArbSeries,AcbSeries})
    # IMPROVE: From the differential equation we could get a
    # recurrence relation for the coefficients. This should be more
    # efficient.

    z₀ = z[0]

    res = zero(z)

    for n = 0:Arblib.degree(z)
        res[n] = U_dz(a, b, z₀, n) / factorial(n)
    end

    return ArbExtras.compose_zero!(res, res, z)
end

U(a::T, b::T, z::T) where {T<:Union{Float64,ComplexF64}} = Arblib.fpwrap_hypgeom_u(a, b, z)

# Used for differentiation w.r.t. a variable which both a and z depend
# on.
function U(a::T, b, z::T) where {T<:Union{ArbSeries,AcbSeries}}
    Arblib.degree(a) == Arblib.degree(z) == 1 ||
        throw(ArgumentError("only supports degree 1"))

    T((U(a[0], b, z[0]), U_da(a[0], b, z[0]) * a[1] + U_dz(a[0], b, z[0]) * z[1]))
end

"""
    U_dz(a, b, z, n::Integer = 1)

Compute ``U(a, b, z)`` differentiated `n` times w.r.t. `z`.

Uses the formula
```
(-1)^n * U(a + n, b + n, z) * rising(a, n)
```
"""
U_dz(a::T, b::T, z::T, n::Integer = 1) where {T} =
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return U(a, b, z)
    else
        return (-1)^n * U(a + n, b + n, z) * rising(a, n)
    end

U_dz(a::Acb, b::Acb, z::AcbSeries, n::Integer = 1) =
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return U(a, b, z)
    else
        return (-1)^n * U(a + n, b + n, z) * rising(a, n)
    end

U_dz(a, b, z, n::Integer = 1) = U_dz(promote(a, b, z)..., n)

"""
    U_da(a, b, z, n::Integer = 1)

Compute ``U(a, b, z)`` differentiated `n` times w.r.t. `a`.

Currently only supports `n <= 1`.

For `n = 1` it uses the asymptotic expansion with a bound for the
remainder term.
"""
function U_da(a::Acb, b::Acb, z::Acb, n::Integer = 1)
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return U(a, b, z)
    elseif n == 1
        abs(z) < 10 && @warn "U_da doesn't work well for small z"

        # IMPROVE: Tune choice of N and improve performance.

        N = 20

        S1 = -sum(0:N-1) do k
            p_U(k, a, b, z)
        end

        S2 = sum(0:N-1) do k
            p_U_da(k, a, b, z)
        end

        R = add_error(
            Acb(0),
            (1 + abs(digamma(a) / log(z))) * C_R_U(N, a, b, z) +
            C_R_U_1(N, a, b, z) +
            C_R_U_2(N, a, b, z),
        )

        return (S1 * log(z) + S2 + R * log(z) * z^-N) * z^-a
    else
        error("no implementation of U_da for n > 1")
    end
end

U_da(a::T, b::T, z::T, n::Integer = 1) where {T} =
    if n == 0
        U(a, b, z)
    else
        res = U_da(Acb(a), Acb(b), Acb(z), n)
        if T <: Real
            Arblib.contains_zero(imag(res)) || error("expected a real result, got $res")
            if T == Float64 && Arblib.rel_accuracy_bits(res) < 50
                @warn "low precision when computing U_da" real(res)
            end
            return convert(T, real(res))
        else
            if T == ComplexF64 && Arblib.rel_accuracy_bits(real(res)) < 50
                @warn "low precision when computing U_da" real(res)
            end
            convert(T, res)
        end
    end
U_da(a, b, z, n::Integer = 1) = U_da(promote(a, b, z)..., n)

"""
    _U_da_finite_difference(a, b, z)

Compute ``U(a, b, z)`` differentiated once w.r.t. `a` using a finite
difference method.

An approximation of the derivative is computed using the central
finite difference
```
(U(a + h, b, z) - U(a - h, b, z)) / 2h
```
for some `h`.

Using Taylor's theorem we get that the error for the approximation is
```
h^2 / 12 * (U_da3(a₁, b, z) + U_da3(a₂, b, z))
```
where we use `U_da3` to denote the third derivative w.r.t. `a`
and `a₁` is a point lying between `a` and `a + h` and `a₂` is a point
lying between `a - h` and `a`.

To bound the third derivative we make use of Cauchy's formula. For `r
> 0` the third derivative at `a` is bounded by
```
6r^(-3) * C
```
where `C` bounds the magnitude of `U(a, b, z)` in the ball of
radius `r` centered around `a`. For the derivatives at `a₁` and `a₂`
we instead get the bound
```
6(r - h)^(-3) * C
```
with for the same `C`.

Combining this we get that the error is bounded by
```
h^2 / (r - h)^3 * C
```
"""
function _U_da_finite_difference(a::Acb, b::Acb, z::Acb)
    # General guidance is to take h to be the square root of the ulp
    # of the argument. In this case we take the ulp to determined by
    # the relative precision of a.
    prec = min(precision(a), Arblib.rel_accuracy_bits(a))
    h = sqrt(eps(Arb(abs(a); prec)))

    # Take a to be somewhat large but so that we avoid the origin
    r = abs(a) / 2

    if !(r > h)
        @warn "Don't have r > h - required for finite difference" r h
        return indeterminate(a)
    end

    res = (U(a + h, b, z) - U(a - h, b, z)) / 2h

    # TODO: Check that this is correct
    C = abs(U(add_error(a, r), b, z))

    error = h^2 / (r - h)^3 * C

    return add_error(res, error)
end

"""
    U_dzda(a, b, z, n::Integer = 1)

Compute ``U(a, b, z)`` differentiated `n` w.r.t. `z` and once w.r.t.
`a`.

For differentiation w.r.t. `z` it uses the formula
```
(-1)^n * U(a + n, b + n, z) * rising(a, n)
```
Which after differentiation gives us
```
(-1)^n * (U(a + n, b + n, z) * drising(a, n) + a * U_da(a + n, b + n, z) * rising(a, n))
```
Where we use `drising(a, n)` to denote the derivative of `rising(a,
n)` w.r.t. `a`.
"""
U_dzda(a::T, b::T, z::T, n::Integer = 1) where {T} =
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return U_da(a, b, z)
    else
        drising = if n == 1 # rising(a, 1) = a
            one(a)
        elseif n == 2 # rising(a, 2) = a * (a + 1)
            2a + 1
        elseif a isa Arb
            rising(ArbSeries((a, 1)), n)[1]
        elseif a isa Acb
            rising(AcbSeries((a, 1)), n)[1]
        else
            throw(ArgumentError("n > 2 only supported for Arb and Acb"))
        end

        res = (-1)^n * (U(a + n, b + n, z) * drising + U_da(a + n, b + n, z) * rising(a, n))

        res
    end

U_dzda(a, b, z, n::Integer = 1) = U_dzda(promote(a, b, z)..., n)

"""
    U_asym_approx(a, b, z)

Compute ``U(a, b, z)``, the confluent hypergeometric function of the
second kind, using the asymptotic series. No attempt is made to bound
the remainder term.
"""
function U_asym_approx(a, b, z)
    N = 20
    res = sum(0:N-1) do k
        rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-z)^k)
    end
    return z^-a * res
end

"""
    U_dz_asym_approx(a, b, z)

Compute ``U(a, b, z)``, the confluent hypergeometric function of the
second kind, differentiated `n` times w.r.t. `z`. Uses the asymptotic
series. No attempt is made to bound the remainder term.

Uses the formula
```
(-1)^n * U(a + n, b + n, z) * rising(a, n)
```
"""
U_dz_asym_approx(a, b, z, n::Integer = 1) =
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return U_asym_approx(a, b, z)
    else
        return (-1)^n * U_asym_approx(a + n, b + n, z) * rising(a, n)
    end
