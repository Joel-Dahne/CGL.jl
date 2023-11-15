"""
    solution_infinity(γ, κ, ξ₁, λ::AbstractGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes `[Q(ξ₁), d(Q)(ξ₁)]`.
"""
function solution_infinity(γ::Acb, κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb}; order = 2)
    v = Arb(0.1) # TODO: How to pick this?

    (; σ, ω, d) = λ

    _, _, c = _abc(κ, λ)

    norm_u = solution_infinity_fixed_point(γ, κ, ξ₁, v, λ)[1]
    norm_u_dξ =
        C_P_dξ(κ, λ, ξ₁) * abs(γ) * ξ₁^(-v - 1) +
        C_u_dξ(κ, ξ₁, v, λ) * norm_u^(2σ + 1) * ξ₁^(2σ * v - 1)
    norm_u_dξ_dξ =
        C_P_dξ_dξ(κ, λ, ξ₁) * abs(γ) * ξ₁^(-v - 2) +
        (
            C_u_dξ_dξ_1(κ, ξ₁, v, λ) * norm_u * ξ₁^(-1) +
            C_u_dξ_dξ_2(κ, ξ₁, v, λ) * norm_u_dξ
        ) *
        norm_u^2σ *
        ξ₁^(2σ * v - 1)

    p = P(ξ₁, (λ, κ))
    p_dξ = P_dξ(ξ₁, (λ, κ))
    e = E(ξ₁, (λ, κ))
    e_dξ = E_dξ(ξ₁, (λ, κ))

    if order == 1
        I_E = I_E_0(γ, κ, ξ₁, v, norm_u, λ)
        I_E_dξ = I_E_dξ_0(γ, κ, ξ₁, v, norm_u, λ)

        I_P = I_P_0(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ)
        I_P_dξ = I_P_dξ_0(γ, κ, ξ₁, v, norm_u, λ)
    elseif order == 2
        I_P_main = let p0 = p_P(0, κ, λ), h = Acb(-2 / σ + d - 3, -2ω / κ)
            abs(γ * p0)^2σ * γ * p0 * B_W(κ, λ) * p0 / 2 *
            ξ₁^(1 + h) *
            expint((1 - h) / 2, c * ξ₁^2)
        end

        I_E_dξ_main = abs(γ)^2σ * γ * J_E(ξ₁, (λ, κ)) * abs(p)^2σ * p
        I_P_dξ_main = -abs(γ)^2σ * γ * J_P(ξ₁, (λ, κ)) * abs(p)^2σ * p

        I_E = zero(γ)
        I_P = I_P_main + I_P_R(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ)

        I_E_dξ = I_E_dξ_main + I_E_dξ_R(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ)
        I_P_dξ = I_P_dξ_main + I_P_dξ_R(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ)
    else
        throw(ArgumentError("invalid value $order for order"))
    end

    Q = γ * p + p * I_E + e * I_P
    dQ = γ * p_dξ + p_dξ * I_E + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ

    return SVector(Q, dQ)
end

function solution_infinity(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::AbstractGLParams{Float64};
    order = 2,
)
    (; σ, ω, d) = λ

    _, _, c = _abc(κ, λ)

    p = P(ξ₁, (λ, κ))
    p_dξ = P_dξ(ξ₁, (λ, κ))
    e = E(ξ₁, (λ, κ))
    e_dξ = E_dξ(ξ₁, (λ, κ))

    if order == 1
        I_E = zero(γ)
        I_P = zero(γ)

        I_E_dξ = zero(γ)
        I_P_dξ = zero(γ)
    elseif order == 2
        I_E = zero(γ)
        I_P = let p0 = p_P(0, κ, λ), h = -2 / σ + d - 3 - im * 2ω / κ
            abs(γ * p0)^2σ * γ * p0 * B_W(κ, λ) * p0 / 2 *
            ξ₁^(1 + h) *
            expint((1 - h) / 2, c * ξ₁^2)
        end

        I_E_dξ = abs(γ)^2σ * γ * J_E(ξ₁, (λ, κ)) * abs(p)^2σ * p
        I_P_dξ = -abs(γ)^2σ * γ * J_P(ξ₁, (λ, κ)) * abs(p)^2σ * p
    else
        throw(ArgumentError("invalid value $order for order"))
    end

    Q = γ * p + p * I_E + e * I_P
    dQ = γ * p_dξ + p_dξ * I_E + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ

    return SVector(Q, dQ)
end

"""
    solution_infinity_asym(γ, κ, ξ₁, λ::AbstractGLParams)

Similar to [`solution_infinity`](@ref) but uses the two leading terms
in the asymptotic expansion directly.
"""
function solution_infinity_asym(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::AbstractGLParams{Float64};
    order = 2,
)
    (; σ, ϵ, δ) = λ

    a, b, c = _abc(κ, λ)
    p0 = p_P(0, κ, λ)

    if order == 1
        γ₁ = 0

        a0 = (γ + γ₁) * p0

        Q = a0 * ξ₁^(-2a)
        dQ = (-2a) * a0 * ξ₁^(-2a - 1)
    elseif order == 2
        γ₁ = abs(γ * p0)^2σ * γ * p0 * (B_W(κ, λ) * p_E(0, κ, λ)) / 2 * ξ₁^-2

        p1 = p_P(1, κ, λ)

        a0 = (γ + γ₁) * p0
        a1 = a0 * (4a * (a - b + 1) * (1 - im * ϵ) + (1 + im * δ) * abs(a0)^2σ) / (2im * κ)

        Q = (a0 + a1 * ξ₁^-2) * ξ₁^(-2a)
        dQ = ((-2a) * a0 + (-2a - 2) * a1 * ξ₁^-2) * ξ₁^(-2a - 1)
    else
        throw(ArgumentError("invalid value $order for order"))
    end

    return SVector(Q, dQ)
end

"""
    solution_infinity_jacobian(γ, κ, ξ₁, λ::AbstractGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes Jacobian w.r.t. the parameters `γ` and `κ` of
`[Q(ξ₁), d(Q)(ξ₁)]. The Jacobian is given by
```
[
d(Q(ξ₁), μ) d(Q(ξ₁), κ)
d(d(Q)(ξ₁), μ) d((Q)(ξ₁), κ)
]
```
where we use `d(Q, μ)` to denote the derivative of `Q` w.r.t. `μ`.
"""
function solution_infinity_jacobian(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    λ::AbstractGLParams{Arb};
    order = 2,
)
    v = Arb(0.1) # TODO: How to pick this?

    (; σ, ω, d) = λ

    _, _, c = _abc(κ, λ)

    norm_u = solution_infinity_fixed_point(γ, κ, ξ₁, v, λ)[1]

    norm_u_dξ =
        C_P_dξ(κ, λ, ξ₁) * abs(γ) * ξ₁^(-v - 1) +
        C_u_dξ(κ, ξ₁, v, λ) * norm_u^(2σ + 1) * ξ₁^(2σ * v - 1)
    norm_u_dξ_dξ =
        C_P_dξ_dξ(κ, λ, ξ₁) * abs(γ) * ξ₁^(-v - 2) +
        (
            C_u_dξ_dξ_1(κ, ξ₁, v, λ) * norm_u * ξ₁^(-1) +
            C_u_dξ_dξ_2(κ, ξ₁, v, λ) * norm_u_dξ
        ) *
        norm_u^2σ *
        ξ₁^(2σ * v - 1)
    norm_u_dγ = let
        CT1, CT2 = C_T1(v, κ, λ, ξ₁)
        num = CT1 * ξ₁^-v
        den = (1 - (2σ + 1) * CT2 * ξ₁^(-2 + 2σ * v) * norm_u^2σ)

        if Arblib.ispositive(den)
            num / den
        else
            @warn "Not positive denominator for norm_u_dγ" num den
            indeterminate(num)
        end
    end
    norm_u_dξ_dγ =
        C_P_dξ(κ, λ, ξ₁) * ξ₁^(-v - 1) +
        (2σ + 1) * C_u_dξ(κ, ξ₁, v, λ) * norm_u^2σ * norm_u_dγ * ξ₁^(2λ.σ * v - 1)
    norm_u_dκ = let
        num = (
            C_u_dκ_1(κ, ξ₁, v, λ) * abs(γ) +
            (
                C_u_dκ_2(κ, ξ₁, v, λ) * norm_u^2 +
                C_u_dκ_3(κ, ξ₁, v, λ) * norm_u * norm_u_dξ +
                C_u_dκ_4(κ, ξ₁, v, λ) * norm_u_dξ^2 +
                C_u_dκ_5(κ, ξ₁, v, λ) * norm_u * norm_u_dξ_dξ
            ) * norm_u^(2σ - 1)
        )
        den = (1 - C_u_dκ_6(κ, ξ₁, v, λ) * norm_u^2σ)

        if Arblib.ispositive(den)
            num / den
        else
            @warn "Not positive denominator for norm_u_dκ" num den
            indeterminate(num)
        end
    end

    if order == 1
        I_E = I_E_0(γ, κ, ξ₁, v, norm_u, λ)
        I_P = I_P_0(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ)
        I_E_dξ = I_E_dξ_0(γ, κ, ξ₁, v, norm_u, λ)
        I_P_dξ = I_P_dξ_0(γ, κ, ξ₁, v, norm_u, λ)

        I_E_dγ = I_E_dγ_0(γ, κ, ξ₁, v, norm_u, norm_u_dγ, λ)
        I_P_dγ = I_P_dγ_0(γ, κ, ξ₁, v, norm_u, norm_u_dγ, norm_u_dξ, norm_u_dξ_dγ, λ)
        I_E_dξ_dγ = I_E_dξ_dγ_0(γ, κ, ξ₁, v, norm_u, norm_u_dγ, λ)
        I_P_dξ_dγ = I_P_dξ_dγ_0(γ, κ, ξ₁, v, norm_u, norm_u_dγ, λ)

        I_E_dκ = I_E_dκ_0(γ, κ, ξ₁, v, norm_u, norm_u_dκ, λ)
        I_P_dκ = I_P_dκ_0(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dκ, λ)
        I_E_dξ_dκ = I_E_dξ_dκ_0(γ, κ, ξ₁, v, norm_u, norm_u_dκ, λ)
        I_P_dξ_dκ = I_P_dξ_dκ_0(γ, κ, ξ₁, v, norm_u, norm_u_dκ, λ)
    elseif order == 2
        I_E_dξ_main, I_P_dξ_main = let p = P(ξ₁, (λ, κ))
            I_E_dξ_main = abs(γ)^2σ * γ * J_E(ξ₁, (λ, κ)) * abs(p)^2σ * p
            I_P_dξ_main = -abs(γ)^2σ * γ * J_P(ξ₁, (λ, κ)) * abs(p)^2σ * p

            I_E_dξ_main, I_P_dξ_main
        end

        I_E = I_E_0(γ, κ, ξ₁, v, norm_u, λ)
        I_P = I_P_0(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ)
        I_E_dξ = I_E_dξ_main + I_E_dξ_R(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ)
        I_P_dξ = I_P_dξ_main + I_P_dξ_R(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ)

        # No higher order bounds for these
        I_E_dγ = I_E_dγ_0(γ, κ, ξ₁, v, norm_u, norm_u_dγ, λ)
        I_P_dγ = I_P_dγ_0(γ, κ, ξ₁, v, norm_u, norm_u_dγ, norm_u_dξ, norm_u_dξ_dγ, λ)
        I_E_dξ_dγ = I_E_dξ_dγ_0(γ, κ, ξ₁, v, norm_u, norm_u_dγ, λ)
        I_P_dξ_dγ = I_P_dξ_dγ_0(γ, κ, ξ₁, v, norm_u, norm_u_dγ, λ)

        # No higher order bounds for these
        I_E_dκ = I_E_dκ_0(γ, κ, ξ₁, v, norm_u, norm_u_dκ, λ)
        I_P_dκ = I_P_dκ_0(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dκ, λ)
        I_E_dξ_dκ = I_E_dξ_dκ_0(γ, κ, ξ₁, v, norm_u, norm_u_dκ, λ)
        I_P_dξ_dκ = I_P_dξ_dκ_0(γ, κ, ξ₁, v, norm_u, norm_u_dκ, λ)
    else
        throw(ArgumentError("invalid value $order for order"))
    end

    Q_dγ = P(ξ₁, (λ, κ)) * (one(γ) + I_E_dγ) + E(ξ₁, (λ, κ)) * I_P_dγ

    dQ_dγ =
        P(ξ₁, (λ, κ)) * I_E_dξ_dγ +
        P_dξ(ξ₁, (λ, κ)) * (one(γ) + I_E_dγ) +
        E(ξ₁, (λ, κ)) * I_P_dξ_dγ +
        E_dξ(ξ₁, (λ, κ)) * I_P_dγ

    Q_dκ =
        P(ξ₁, (λ, κ)) * I_E_dκ +
        P_dκ(ξ₁, (λ, κ)) * (γ + I_E) +
        E(ξ₁, (λ, κ)) * I_P_dκ +
        E_dκ(ξ₁, (λ, κ)) * I_P

    dQ_dκ =
        P(ξ₁, (λ, κ)) * I_E_dξ_dκ +
        P_dξ(ξ₁, (λ, κ)) * I_E_dκ +
        P_dκ(ξ₁, (λ, κ)) * I_E_dξ +
        P_dξ_dκ(ξ₁, (λ, κ)) * (γ + I_E) +
        E(ξ₁, (λ, κ)) * I_P_dξ_dκ +
        E_dξ(ξ₁, (λ, κ)) * I_P_dκ +
        E_dκ(ξ₁, (λ, κ)) * I_P_dξ +
        E_dξ_dκ(ξ₁, (λ, κ)) * I_P

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end

function solution_infinity_jacobian(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::AbstractGLParams{Float64},
)
    # FIXME: Take into account the rest of Q
    Q_dγ = P(ξ₁, (λ, κ))
    dQ_dγ = P_dξ(ξ₁, (λ, κ))

    Q_dκ = γ * P_dκ(ξ₁, (λ, κ))
    dQ_dκ = γ * P_dξ_dκ(ξ₁, (λ, κ))

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end
