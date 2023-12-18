"""
    solution_infinity(γ, κ, ξ₁, λ::CGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes `[Q(ξ₁), d(Q)(ξ₁)]`.
"""
function solution_infinity(γ::Acb, κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    v = Arb(0.1) # TODO: How to pick this?

    (; σ) = λ

    norm_u = norm_bound_u(γ, κ, ξ₁, v, λ)
    norm_u_dξ = norm_bound_u_dξ(γ, κ, ξ₁, v, norm_u, λ)
    norm_u_dξ_dξ = norm_bound_u_dξ_dξ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ)
    norm_u_dξ_dξ_dξ = norm_bound_u_dξ_dξ_dξ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ)

    # Compute zeroth order bounds
    Q = add_error(zero(γ), norm_u * ξ₁^(-1 / σ + v))
    dQ = add_error(zero(γ), norm_u_dξ * ξ₁^(-1 / σ + v))

    # Improve the bounds iteratively

    p = P(ξ₁, (λ, κ))
    p_dξ = P_dξ(ξ₁, (λ, κ))
    e = E(ξ₁, (λ, κ))
    e_dξ = E_dξ(ξ₁, (λ, κ))

    # IMPROVE: In practice three iterations seems to be enough to
    # saturate the convergence. But it might be better to choose this
    # dynamically.
    for _ = 1:3
        I_E = zero(γ)
        I_P = I_P_1(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dξ_dξ_dξ, Q, dQ, λ)

        Q = γ * p + p * I_E + e * I_P

        I_E_dξ = J_E(ξ₁, (λ, κ)) * abs(Q)^2σ * Q
        I_P_dξ = -J_P(ξ₁, (λ, κ)) * abs(Q)^2σ * Q

        dQ = γ * p_dξ + p_dξ * I_E + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ
    end

    return SVector(Q, dQ)
end

function solution_infinity(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64};
    order = 2,
)
    (; σ, ω, d) = λ

    _, _, c = _abc(κ, λ)

    p = P(ξ₁, (λ, κ))
    p_dξ = P_dξ(ξ₁, (λ, κ))
    e = E(ξ₁, (λ, κ))
    e_dξ = E_dξ(ξ₁, (λ, κ))

    # Compute first order approximation
    Q_1 = γ * p
    dQ_1 = γ * p_dξ

    order == 1 && return SVector(Q_1, dQ_1)

    # Compute second order approximation
    Q_2, dQ_2 = let
        I_E = zero(γ)
        I_P = let p0 = p_P(0, κ, λ), h = -2 / σ + d - 3 - im * 2ω / κ
            abs(γ * p0)^2σ * γ * p0 * B_W(κ, λ) * p0 / 2 *
            ξ₁^(1 + h) *
            expint((1 - h) / 2, c * ξ₁^2)
        end

        Q = γ * p + p * I_E + e * I_P

        I_E_dξ = J_E(ξ₁, (λ, κ)) * abs(Q)^2σ * Q
        I_P_dξ = -J_P(ξ₁, (λ, κ)) * abs(Q)^2σ * Q

        dQ = γ * p_dξ + p_dξ * I_E + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ

        Q, dQ
    end

    order == 2 && return SVector(Q_2, dQ_2)

    throw(ArgumentError("invalid value $order for order"))
end

"""
    solution_infinity_asym(γ, κ, ξ₁, λ::CGLParams)

Similar to [`solution_infinity`](@ref) but uses the two leading terms
in the asymptotic expansion directly.
"""
function solution_infinity_asym(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64};
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
    solution_infinity_jacobian(γ, κ, ξ₁, λ::CGLParams)

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
function solution_infinity_jacobian(γ::Acb, κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    v = Arb(0.1) # TODO: How to pick this?

    (; σ) = λ

    norm_u = norm_bound_u(γ, κ, ξ₁, v, λ)
    norm_u_dξ = norm_bound_u_dξ(γ, κ, ξ₁, v, norm_u, λ)
    norm_u_dξ_dξ = norm_bound_u_dξ_dξ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ)
    norm_u_dξ_dξ_dξ = norm_bound_u_dξ_dξ_dξ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ)

    norm_u_dγ = norm_bound_u_dγ(γ, κ, ξ₁, v, norm_u, λ)
    norm_u_dξ_dγ = norm_bound_u_dξ_dγ(γ, κ, ξ₁, v, norm_u, norm_u_dγ, λ)
    norm_u_dκ = norm_bound_u_dκ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ)
    norm_u_dξ_dκ =
        norm_bound_u_dξ_dκ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dκ, λ)

    # Compute zeroth order bounds
    Q = add_error(zero(γ), norm_u * ξ₁^(-1 / σ + v))
    dQ = add_error(zero(γ), norm_u_dξ * ξ₁^(-1 / σ + v))
    Q_dγ = add_error(zero(γ), norm_u_dγ * ξ₁^(-1 / σ + v))
    dQ_dγ = add_error(zero(γ), norm_u_dξ_dγ * ξ₁^(-1 / σ + v))
    Q_dκ = add_error(zero(γ), norm_u_dκ * ξ₁^(-1 / σ + v))
    dQ_dκ = add_error(zero(γ), norm_u_dξ_dκ * ξ₁^(-1 / σ + v))

    # Improve the bounds iteratively

    p = P(ξ₁, (λ, κ))
    p_dξ = P_dξ(ξ₁, (λ, κ))
    p_dκ = P_dκ(ξ₁, (λ, κ))
    p_dξ_dκ = P_dξ_dκ(ξ₁, (λ, κ))
    e = E(ξ₁, (λ, κ))
    e_dξ = E_dξ(ξ₁, (λ, κ))
    e_dκ = E_dκ(ξ₁, (λ, κ))
    e_dξ_dκ = E_dξ_dκ(ξ₁, (λ, κ))

    # IMPROVE: In practice three iterations seems to be enough to
    # saturate the convergence. But it might be better to choose this
    # dynamically.
    for _ = 1:3
        I_E = zero(γ)
        I_P = I_P_1(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dξ_dξ_dξ, Q, dQ, λ)

        I_E_dγ = zero(γ)
        I_P_dγ = I_P_dγ_0(γ, κ, ξ₁, v, norm_u, norm_u_dγ, norm_u_dξ, norm_u_dξ_dγ, λ)

        I_E_dκ = zero(γ)
        I_P_dκ = I_P_dκ_0(
            γ,
            κ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dξ_dξ,
            norm_u_dκ,
            norm_u_dξ_dκ,
            λ,
        )

        Q = γ * p + p * I_E + e * I_P
        Q_dγ = p + p * I_E_dγ + e * I_P_dγ
        Q_dκ = γ * p_dκ + p * I_E_dκ + p_dκ * I_E + e * I_P_dκ + e_dκ * I_P

        I_E_dξ = J_E(ξ₁, (λ, κ)) * abs(Q)^2σ * Q
        I_P_dξ = -J_P(ξ₁, (λ, κ)) * abs(Q)^2σ * Q

        I_E_dξ_dγ =
            J_E(ξ₁, (λ, κ)) *
            abs(Q)^(2σ - 2) *
            (2σ * real(conj(Q) * Q_dγ) * Q + abs(Q)^2 * Q_dγ)
        I_P_dξ_dγ =
            -J_P(ξ₁, (λ, κ)) *
            abs(Q)^(2σ - 2) *
            (2σ * real(conj(Q) * Q_dγ) * Q + abs(Q)^2 * Q_dγ)

        I_E_dξ_dκ =
            J_E_dκ(ξ₁, (λ, κ)) * abs(Q)^2σ * Q +
            J_E(ξ₁, (λ, κ)) *
            abs(Q)^(2σ - 2) *
            (2σ * real(conj(Q) * Q_dκ) * Q + abs(Q)^2 * Q_dκ)
        I_P_dξ_dκ =
            -J_P_dκ(ξ₁, (λ, κ)) * abs(Q)^2σ * Q -
            J_P(ξ₁, (λ, κ)) *
            abs(Q)^(2σ - 2) *
            (2σ * real(conj(Q) * Q_dκ) * Q + abs(Q)^2 * Q_dκ)

        dQ = γ * p_dξ + p_dξ * I_E + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ
        dQ_dγ = p_dξ + p * I_E_dξ_dγ + p_dξ * I_E_dγ + e * I_P_dξ_dγ + e_dξ * I_P_dγ
        dQ_dκ =
            γ * p_dξ_dκ +
            p * I_E_dξ_dκ +
            p_dξ * I_E_dκ +
            p_dκ * I_E_dξ +
            p_dξ_dκ * I_E +
            e * I_P_dξ_dκ +
            e_dξ * I_P_dκ +
            e_dκ * I_P_dξ +
            e_dξ_dκ * I_P
    end

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end

function solution_infinity_jacobian(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    # IMPROVE: Add higher order versions
    Q_dγ = P(ξ₁, (λ, κ))
    dQ_dγ = P_dξ(ξ₁, (λ, κ))

    Q_dκ = γ * P_dκ(ξ₁, (λ, κ))
    dQ_dκ = γ * P_dξ_dκ(ξ₁, (λ, κ))

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end
