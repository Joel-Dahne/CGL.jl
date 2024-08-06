"""
    cgl_equation_real(Q, κ, ϵ, ξ, λ)
    cgl_equation_real(Q, (κ, ϵ, λ), ξ)

Evaluate the right hand side of the ODE when written as a four
dimensional real system. It is evaluated at the point
```
Q = [a, b, α, β]
```
and time `ξ`.

For `λ.d != 1` there is a removable singularity at `ξ = 0`. To return
a finite value we in this case required that `α = β = 0`.
"""
function cgl_equation_real(Q, κ, ϵ, ξ, λ::CGLParams)
    (; d, ω, σ, δ) = λ
    a, b, α, β = Q

    a2b2σ = (a^2 + b^2)^σ

    # TODO: Should we keep fastmath? It makes it faster and probably
    # improves accuracy due to fma.
    @fastmath F1 = κ * ξ * β + κ / σ * b + ω * a - a2b2σ * a + δ * a2b2σ * b
    @fastmath F2 = -κ * ξ * α - κ / σ * a + ω * b - a2b2σ * b - δ * a2b2σ * a

    if !isone(d) && !(iszero(ξ) && iszero(α) && iszero(β))
        F1 -= (d - 1) / ξ * (α + ϵ * β)
        F2 -= (d - 1) / ξ * (β - ϵ * α)
    end

    return SVector(α, β, (F1 - ϵ * F2) / (1 + ϵ^2), (ϵ * F1 + F2) / (1 + ϵ^2))
end

# For use with ODEProblem
cgl_equation_real(u, (κ, ϵ, λ), ξ) = cgl_equation_real(u, κ, ϵ, ξ, λ)

"""
    cgl_equation_real_taylor(((a0, a1), (b0, b1)), κ, ϵ, ξ₀, λ; degree = 5)

Compute the Taylor expansions of `a` and `b` in
[`cgl_equation_real`](@ref). The expansions are centered at the point
`ξ = ξ₀`, with the first two coefficients in the expansions for `a`
and `b` given by `a0, a1` and `b0, b1` respectively.
"""
function cgl_equation_real_taylor(
    Q_ξ₀::AbstractVector{NTuple{2,Arb}},
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree::Integer = 5,
)
    (; d, ω, σ, δ) = λ

    a = ArbSeries(Q_ξ₀[1]; degree)
    b = ArbSeries(Q_ξ₀[2]; degree)

    if iszero(ξ₀) && !isone(d) && !iszero(a[1]) && !iszero(b[1])
        return SVector(indeterminate(a), indeterminate(b))
    end

    for n = 0:degree-2
        a2b2σ = (ArbSeries(a, degree = n + 1)^2 + ArbSeries(b, degree = n + 1)^2)^σ
        u1 = a2b2σ * a
        u2 = a2b2σ * b

        if iszero(ξ₀)
            F1 = κ * n * b[n] + κ / σ * b[n] + ω * a[n] - u1[n] + δ * u2[n]
            F2 = -κ * n * a[n] - κ / σ * a[n] + ω * b[n] - u2[n] - δ * u1[n]

            a[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + d) * (1 + ϵ^2))
            b[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + d) * (1 + ϵ^2))
        else
            F1 =
                κ * (ξ₀ * (n + 1) * b[n+1] + n * b[n]) + κ / σ * b[n] + ω * a[n] - u1[n] +
                δ * u2[n]
            F2 =
                -κ * (ξ₀ * (n + 1) * a[n+1] + n * a[n]) - κ / σ * a[n] + ω * b[n] - u2[n] -
                δ * u1[n]

            if !isone(d)
                v1 = Arblib.derivative(a) / ArbSeries((ξ₀, 1), degree = n) # a' / ξ
                v2 = Arblib.derivative(b) / ArbSeries((ξ₀, 1), degree = n) # b' / ξ
                F1 += (1 - d) * (v1[n] + ϵ * v2[n])
                F2 += (1 - d) * (v2[n] - ϵ * v1[n])
            end

            a[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + 1) * (1 + ϵ^2))
            b[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + 1) * (1 + ϵ^2))
        end
    end

    return SVector(a, b)
end

"""
    cgl_equation_real_dμ_taylor(((a0_dμ, a1_dμ), (b0_dμ, b1_dμ)), a, b, κ, ϵ, ξ₀, λ; degree = 5)

Compute the Taylor expansions of `a` and `b` in
[`cgl_equation_real`](@ref) differentiated with respect to `μ`. The
expansions are centered at the point `ξ = ξ₀` with the first two
coefficients in the expansions given by `a0_dμ, a1_dμ` and `b0_dμ,
b1_dμ` respectively.
"""
function cgl_equation_real_dμ_taylor(
    Q_ξ₀_dμ::AbstractVector{NTuple{2,Arb}},
    a::ArbSeries,
    b::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree::Integer = 5,
)
    @assert Arblib.degree(a) == Arblib.degree(b) == degree

    (; d, ω, σ, δ) = λ

    a_dμ = ArbSeries(Q_ξ₀_dμ[1]; degree)
    b_dμ = ArbSeries(Q_ξ₀_dμ[2]; degree)

    if iszero(ξ₀) && !isone(d) && !iszero(a_dμ[1]) && !iszero(b_dμ[1])
        return SVector(indeterminate(a_dμ), indeterminate(b_dμ))
    end

    a2b2σ = (a^2 + b^2)^σ
    d_a2b2σ_aa = 2σ * (a^2 + b^2)^(σ - 1) * a^2
    d_a2b2σ_ab = 2σ * (a^2 + b^2)^(σ - 1) * a * b
    d_a2b2σ_bb = 2σ * (a^2 + b^2)^(σ - 1) * b^2

    for n = 0:degree-2
        u1 =
            a2b2σ * ArbSeries(a_dμ, degree = n + 1) +
            d_a2b2σ_aa * ArbSeries(a_dμ, degree = n + 1) +
            d_a2b2σ_ab * ArbSeries(b_dμ, degree = n + 1)
        u2 =
            a2b2σ * ArbSeries(b_dμ, degree = n + 1) +
            d_a2b2σ_ab * ArbSeries(a_dμ, degree = n + 1) +
            d_a2b2σ_bb * ArbSeries(b_dμ, degree = n + 1)

        if iszero(ξ₀)
            F1 = κ * n * b_dμ[n] + κ / σ * b_dμ[n] + ω * a_dμ[n] - u1[n] + δ * u2[n]
            F2 = -κ * n * a_dμ[n] - κ / σ * a_dμ[n] + ω * b_dμ[n] - u2[n] - δ * u1[n]

            a_dμ[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + d) * (1 + ϵ^2))
            b_dμ[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + d) * (1 + ϵ^2))
        else
            error("not implemented")
        end
    end

    return SVector(a_dμ, b_dμ)
end

"""
    cgl_equation_real_dκ_taylor(((a0_dμ, a1_dμ), (b0_dμ, b1_dμ)), a, b, κ, ϵ, ξ₀, λ; degree = 5)

Compute the Taylor expansions of `a` and `b` in
[`cgl_equation_real`](@ref) differentiated with respect to `κ`. The
expansions are centered at the point `ξ = ξ₀` with the first two
coefficients in the expansions given by `a0_dκ, a1_dκ` and `b0_dκ,
b1_dκ` respectively.
"""
function cgl_equation_real_dκ_taylor(
    Q_ξ₀_dκ::AbstractVector{NTuple{2,Arb}},
    a::ArbSeries,
    b::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree::Integer = 5,
)
    @assert Arblib.degree(a) == Arblib.degree(b) == degree

    (; d, ω, σ, δ) = λ

    a_dκ = ArbSeries(Q_ξ₀_dκ[1]; degree)
    b_dκ = ArbSeries(Q_ξ₀_dκ[2]; degree)

    if iszero(ξ₀) && !isone(d) && !iszero(a_dκ[1]) && !iszero(b_dκ[1])
        return SVector(indeterminate(a_dκ), indeterminate(b_dκ))
    end

    a2b2σ = (a^2 + b^2)^σ
    d_a2b2σ_aa = 2σ * (a^2 + b^2)^(σ - 1) * a^2
    d_a2b2σ_ab = 2σ * (a^2 + b^2)^(σ - 1) * a * b
    d_a2b2σ_bb = 2σ * (a^2 + b^2)^(σ - 1) * b^2

    for n = 0:degree-2
        u1 =
            a2b2σ * ArbSeries(a_dκ, degree = n + 1) +
            d_a2b2σ_aa * ArbSeries(a_dκ, degree = n + 1) +
            d_a2b2σ_ab * ArbSeries(b_dκ, degree = n + 1)
        u2 =
            a2b2σ * ArbSeries(b_dκ, degree = n + 1) +
            d_a2b2σ_ab * ArbSeries(a_dκ, degree = n + 1) +
            d_a2b2σ_bb * ArbSeries(b_dκ, degree = n + 1)

        if iszero(ξ₀)
            F1 =
                κ * n * b_dκ[n] + n * b[n] + κ / σ * b_dκ[n] + 1 / σ * b[n] + ω * a_dκ[n] -
                u1[n] + δ * u2[n]
            F2 =
                -κ * n * a_dκ[n] - n * a[n] - κ / σ * a_dκ[n] - 1 / σ * a[n] + ω * b_dκ[n] -
                u2[n] - δ * u1[n]

            a_dκ[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + d) * (1 + ϵ^2))
            b_dκ[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + d) * (1 + ϵ^2))
        else
            error("not implemented")
        end
    end

    return SVector(a_dκ, b_dκ)
end

"""
    cgl_equation_real_dϵ_taylor(((a_dϵ0, a_dϵ1), (b_dϵ0, b_dϵ1)), a, b, κ, ϵ, ξ₀, λ; degree = 5)

Compute the Taylor expansions of `a` and `b` in
[`cgl_equation_real`](@ref) differentiated with respect to `ϵ`. The
expansions are centered at the point `ξ = ξ₀` with the first two
coefficients in the expansions given by `a0_dϵ, a1_dϵ` and `b0_dϵ,
b1_dϵ` respectively.
"""
function cgl_equation_real_dϵ_taylor(
    Q_ξ₀_dϵ::AbstractVector{NTuple{2,Arb}},
    a::ArbSeries,
    b::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree::Integer = 5,
)
    @assert Arblib.degree(a) == Arblib.degree(b) == degree

    (; d, ω, σ, δ) = λ

    a_dϵ = ArbSeries(Q_ξ₀_dϵ[1]; degree)
    b_dϵ = ArbSeries(Q_ξ₀_dϵ[2]; degree)

    if iszero(ξ₀) && !isone(d) && !iszero(a_dϵ[1]) && !iszero(b_dϵ[1])
        return SVector(indeterminate(a_dϵ), indeterminate(b_dϵ))
    end

    a2b2σ = (a^2 + b^2)^σ
    d_a2b2σ_aa = 2σ * (a^2 + b^2)^(σ - 1) * a^2
    d_a2b2σ_ab = 2σ * (a^2 + b^2)^(σ - 1) * a * b
    d_a2b2σ_bb = 2σ * (a^2 + b^2)^(σ - 1) * b^2

    for n = 0:degree-2
        u1 =
            a2b2σ * ArbSeries(a_dϵ, degree = n + 1) +
            d_a2b2σ_aa * ArbSeries(a_dϵ, degree = n + 1) +
            d_a2b2σ_ab * ArbSeries(b_dϵ, degree = n + 1)
        u2 =
            a2b2σ * ArbSeries(b_dϵ, degree = n + 1) +
            d_a2b2σ_ab * ArbSeries(a_dϵ, degree = n + 1) +
            d_a2b2σ_bb * ArbSeries(b_dϵ, degree = n + 1)

        if iszero(ξ₀)
            F1 =
                -(n + 2) * (n + d) * b[n+2] +
                κ * n * b_dϵ[n] +
                κ / σ * b_dϵ[n] +
                ω * a_dϵ[n] - u1[n] + δ * u2[n]
            F2 =
                (n + 2) * (n + d) * a[n+2] - κ * n * a_dϵ[n] - κ / σ * a_dϵ[n] +
                ω * b_dϵ[n] - u2[n] - δ * u1[n]

            a_dϵ[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + d) * (1 + ϵ^2))
            b_dϵ[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + d) * (1 + ϵ^2))
        else
            error("not implemented")
        end
    end

    return SVector(a_dϵ, b_dϵ)
end
