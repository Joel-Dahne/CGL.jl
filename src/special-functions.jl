export rising, hypgeom_u, hypgeom_u_dz

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
