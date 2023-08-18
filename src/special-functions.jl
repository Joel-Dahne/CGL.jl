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
rising(x, n::Integer) = prod(x + i for i = 0:n-1)

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
    else
        a = AcbSeries((a, 1), degree = n)
        b = AcbSeries(b, degree = n)
        z = AcbSeries(z, degree = n)
        res = zero(a)
        Arblib.hypgeom_u_1f1_series!(res, a, b, z, n + 1)
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
            return convert(T, real(res))
        else
            convert(T, res)
        end
    end
hypgeom_u_da(a, b, z, n::Integer = 1) = hypgeom_u_da(promote(a, b, z)..., n)

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
