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

U(a::T, b::T, z::T) where {T<:Union{Float64,ComplexF64}} = Arblib.fpwrap_hypgeom_u(a, b, z)

# IMPROVE: From the differential equation we could get a recurrence
# relation for the coefficients. This should be more efficient.
function U(a, b, z::Union{ArbSeries,AcbSeries})
    z₀ = z[0]
    res = zero(z)
    for n = 0:Arblib.degree(z)
        res[n] = U_dz(a, b, z₀, n) / factorial(n)
    end
    return ArbExtras.compose_zero!(res, res, z)
end

# This is used for differentiation w.r.t. a variable which both a and
# z depend on.
function U(a::T, b, z::T) where {T<:Union{ArbSeries,AcbSeries}}
    Arblib.degree(a) == Arblib.degree(z) == 1 ||
        throw(ArgumentError("only supports degree 1"))

    # Constant term is given by U(a[0], b, z[0]).
    # Derivative is given by combination of derivative with respect to
    # a and z.
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
U_dz(a, b, z, n::Integer = 1) = U_dz(promote(a, b, z)..., n)

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

"""
    U_da(a, b, z)

Compute ``U(a, b, z)`` differentiated w.r.t. `a`.

The computation is based on Lemma REF(lemma:U-a), using an asymptotic
expansion combined with a bound for the remainder.
"""
U_da(a, b, z) = U_da(promote(a, b, z)...)

function U_da(a::Acb, b::Acb, z::Acb)
    abs(z) < 10 && @debug "U_da doesn't work well for small z"

    # IMPROVE: Tune choice of N
    N = 20

    # Compute the two sums in the asymptotic expansion
    S1 = zero(a)
    S2 = zero(a)
    term1 = zero(a)
    term2 = zero(a)
    tmp1 = zero(a)
    tmp2 = zero(a)
    tmp3 = zero(a)
    for k = 0:(N-1)
        p_U_p_U_da!(term1, term2, k, a, b, z, tmp1, tmp2, tmp3)
        Arblib.add!(S1, S1, term1)
        Arblib.add!(S2, S2, term2)
    end

    # Compute bounds for the remainder terms. This is based on Lemma
    # REF(lemma:U-a).
    R_U = add_error(Acb(0), C_R_U(N, a, b, z))
    R_U_1, R_U_2 = add_error.(Acb(0), C_R_U_12(N, a, b, z))

    # Compute the final enclosure. Note that we move the
    # multiplication with log(z) inside to cancel some factors.
    # Note that Γ'(a) / Γ(a) is exactly the digamma function.
    return (
        -S1 * log(z) +
        S2 +
        ((1 - digamma(a) / log(z)) * R_U + R_U_1 + R_U_2) * log(z) * z^-N
    ) * z^-a
end

# Converts to Acb and use asymptotic expansion
function U_da(a::T, b::T, z::T) where {T}
    res = U_da(Acb(a), Acb(b), Acb(z))
    if T <: Real
        Arblib.contains_zero(imag(res)) || error("expected a real result, got $res")
        if T == Float64 && Arblib.rel_accuracy_bits(res) < 50
            @debug "low precision when computing U_da" real(res)
        end
        return convert(T, real(res))
    else
        if T == ComplexF64 && Arblib.rel_accuracy_bits(real(res)) < 50
            @debug "low precision when computing U_da" real(res)
        end
        return convert(T, res)
    end
end

# IMPROVE: From the differential equation we could get a recurrence
# relation for the coefficients. This should be more efficient.
function U_da(a, b, z::Union{ArbSeries,AcbSeries})
    z₀ = z[0]
    res = zero(z)
    for n = 0:Arblib.degree(z)
        res[n] = U_dzda(a, b, z₀, n) / factorial(n)
    end
    return ArbExtras.compose_zero!(res, res, z)
end

"""
    U_dzda(a, b, z, n::Integer = 1)

Compute ``U(a, b, z)`` differentiated `n` w.r.t. `z` and once w.r.t.
`a`.

For differentiation w.r.t. `z` it uses the formula
```
(-1)^n * U(a + n, b + n, z) * rising(a, n)
```
Which after differentiation w.r.t. `a` gives us
```
(-1)^n * (U(a + n, b + n, z) * drising(a, n) + U_da(a + n, b + n, z) * rising(a, n))
```
Where we use `drising(a, n)` to denote the derivative of `rising(a,
n)` w.r.t. `a`.
"""
U_dzda(a, b, z, n::Integer = 1) = U_dzda(promote(a, b, z)..., n)

function U_dzda(a::T, b::T, z::T, n::Integer = 1) where {T}
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

        return (-1)^n *
               (U(a + n, b + n, z) * drising + U_da(a + n, b + n, z) * rising(a, n))
    end
end
