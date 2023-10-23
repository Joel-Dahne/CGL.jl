export rising, hypgeom_u, hypgeom_u_dz, hypgeom_u_da, hypgeom_u_dzda

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
    hypgeom_u(a, b, z)

Compute ``U(a, b, z)``, the confluent hypergeometric function of
the second kind.
"""
hypgeom_u(a, b, z) = hypgeom_u(promote(a, b, z)...)

hypgeom_u(a::Arblib.ArbOrRef, b::Arblib.ArbOrRef, z::Arblib.ArbOrRef) =
    Arblib.hypgeom_u!(zero(z), a, b, z)

hypgeom_u(a::Arblib.AcbOrRef, b::Arblib.AcbOrRef, z::Arblib.AcbOrRef) =
    Arblib.hypgeom_u!(zero(z), a, b, z)

function hypgeom_u(a, b, z::Union{ArbSeries,AcbSeries})
    # IMPROVE: From the differential equation we could get a
    # recurrence relation for the coefficients. This should be more
    # efficient.

    z₀ = z[0]

    res = zero(z)

    for n = 0:Arblib.degree(z)
        res[n] = hypgeom_u_dz(a, b, z₀, n) / factorial(n)
    end

    return ArbExtras.compose_zero!(res, res, z)
end

hypgeom_u(a::T, b::T, z::T) where {T<:Union{Float64,ComplexF64}} =
    Arblib.fpwrap_hypgeom_u(a, b, z)

# Used for differentiation w.r.t. a variable which both a and z depend
# on.
function hypgeom_u(a::T, b::Union{Arb,Acb}, z::T) where {T<:Union{ArbSeries,AcbSeries}}
    Arblib.degree(a) == Arblib.degree(z) == 1 ||
        throw(ArgumentError("only supports degree 1"))

    T((
        hypgeom_u(a[0], b, z[0]),
        hypgeom_u_da(a[0], b, z[0]) * a[1] + hypgeom_u_dz(a[0], b, z[0]) * z[1],
    ))
end

"""
    hypgeom_u_dz(a, b, z, n::Integer = 1)

Compute ``U(a, b, z)`` differentiated `n` times w.r.t. `z`.

Uses the formula
```
(-1)^n * hypgeom_u(a + n, b + n, z) * rising(a, n)
```
"""
hypgeom_u_dz(a::T, b::T, z::T, n::Integer = 1) where {T} =
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return hypgeom_u(a, b, z)
    else
        return (-1)^n * hypgeom_u(a + n, b + n, z) * rising(a, n)
    end

hypgeom_u_dz(a::Acb, b::Acb, z::AcbSeries, n::Integer = 1) =
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return hypgeom_u(a, b, z)
    else
        return (-1)^n * hypgeom_u(a + n, b + n, z) * rising(a, n)
    end

hypgeom_u_dz(a, b, z, n::Integer = 1) = hypgeom_u_dz(promote(a, b, z)..., n)

"""
    hypgeom_u_da(a, b, z, n::Integer = 1)

Compute ``U(a, b, z)`` differentiated `n` times w.r.t. `a`.

The computation of the derivative is done using
[`Arblib.hypgeom_u_1f1_series!`](@ref).
"""
function hypgeom_u_da(a::Acb, b::Acb, z::Acb, n::Integer = 1)
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return hypgeom_u(a, b, z)
    elseif n == 1
        return _hypgeom_u_da_finite_difference(a, b, z)
    else
        # TODO: The below code is hopefully not used.

        a_s = AcbSeries((a, 1), degree = n)

        if Arblib.isexact(a) && Arblib.isexact(b) && Arblib.isexact(z)
            b_s = AcbSeries(b, degree = n)
            z_s = AcbSeries(z, degree = n)
            res = zero(a_s)

            # This version only works well at very high precision and
            # not for wide input.
            Arblib.hypgeom_u_1f1_series!(res, a_s, b_s, z_s, n + 1, 10precision(a_s))
        else
            # FIXME: This uses the asymptotic expansion but currently
            # doesn't properly bound the remainder term.
            N = 10
            S = sum(0:N-1) do k
                rising(a_s, k) * rising(a_s - b + 1, k) / (factorial(k) * (-z)^k)
            end
            # Estimated upper bound of error, 10 times the magnitude
            # of the next term.
            term = (rising(a_s, N) * rising(a_s - b + 1, N)) / (factorial(N) * (-z)^N)
            for i = 0:n
                S[i] = add_error(S[i], 10abs(term[i]))
            end

            res = z^-a_s * S
        end

        if n == 1
            return res[1]
        else
            return res[n] * factorial(n)
        end
    end
end

hypgeom_u_da(a::T, b::T, z::T, n::Integer = 1) where {T} =
    if n == 0
        hypgeom_u(a, b, z)
    else
        res = hypgeom_u_da(Acb(a), Acb(b), Acb(z), n)
        if T <: Real
            isreal(res) || error("expected a real result, got $res")
            if T == Float64 && Arblib.rel_accuracy_bits(real(res)) < 50
                @warn "low precision when computing hypgeom_u_da" real(res)
            end
            return convert(T, real(res))
        else
            if T == ComplexF64 && Arblib.rel_accuracy_bits(real(res)) < 50
                @warn "low precision when computing hypgeom_u_da" real(res)
            end
            convert(T, res)
        end
    end
hypgeom_u_da(a, b, z, n::Integer = 1) = hypgeom_u_da(promote(a, b, z)..., n)

"""
    _hypgeom_u_da_finite_difference(a, b, z)

Compute ``U(a, b, z)`` differentiated once w.r.t. `a` using a finite
difference method.

An approximation of the derivative is computed using the central
finite difference
```
(hypgeom_u(a + h, b, z) - hypgeom_u(a - h, b, z)) / 2h
```
for some `h`.

Using Taylor's theorem we get that the error for the approximation is
```
h^2 / 12 * (hypgeom_u_da3(a₁, b, z) + hypgeom_u_da3(a₂, b, z))
```
where we use `hypgeom_u_da3` to denote the third derivative w.r.t. `a`
and `a₁` is a point lying between `a` and `a + h` and `a₂` is a point
lying between `a - h` and `a`.

To bound the third derivative we make use of Cauchy's formula. For `r
> 0` the third derivative at `a` is bounded by
```
6r^(-3) * C
```
where `C` bounds the magnitude of `hypgeom_u(a, b, z)` in the ball of
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
function _hypgeom_u_da_finite_difference(a::Acb, b::Acb, z::Acb)
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

    res = (hypgeom_u(a + h, b, z) - hypgeom_u(a - h, b, z)) / 2h

    # TODO: Check that this is correct
    C = abs(hypgeom_u(add_error(a, r), b, z))

    error = h^2 / (r - h)^3 * C

    return add_error(res, error)
end

"""
    hypgeom_u_dzda(a, b, z)

Compute ``U(a, b, z)`` differentiated once w.r.t. `z` and once w.r.t.
`a`.

For differentiation w.r.t. `z` it uses that the derivative is given by
```
-a * hypgeom_u(a + 1, b + 1, z)
```
Which after differentiation gives us
```
-hypgeom_u(a + 1, b + 1, z) - a * hypgeom_u_da(a + 1, b + 1, z)
```
"""
hypgeom_u_dzda(a::T, b::T, z::T) where {T} =
    -hypgeom_u(a + 1, b + 1, z) - a * hypgeom_u_da(a + 1, b + 1, z)

hypgeom_u_dzda(a, b, z) = hypgeom_u_dzda(promote(a, b, z)...)

"""
    hypgeom_u_asym_approx(a, b, z)

Compute ``U(a, b, z)``, the confluent hypergeometric function of the
second kind, using the asymptotic series. No attempt is made to bound
the remainder term.
"""
function hypgeom_u_asym_approx(a, b, z)
    N = 20
    res = sum(0:N-1) do k
        rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-z)^k)
    end
    return z^-a * res
end

"""
    hypgeom_u_dz_asym_approx(a, b, z)

Compute ``U(a, b, z)``, the confluent hypergeometric function of the
second kind, differentiated `n` times w.r.t. `z`. Uses the asymptotic
series. No attempt is made to bound the remainder term.

Uses the formula
```
(-1)^n * hypgeom_u(a + n, b + n, z) * rising(a, n)
```
"""
hypgeom_u_dz_asym_approx(a, b, z, n::Integer = 1) =
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return hypgeom_u_asym_approx(a, b, z)
    else
        return (-1)^n * hypgeom_u_asym_approx(a + n, b + n, z) * rising(a, n)
    end
