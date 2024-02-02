"""
    solution_infinity(γ, κ, ξ₁, λ::CGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes `[Q(ξ₁), d(Q)(ξ₁)]`.
"""
function solution_infinity(γ::Acb, κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    v = Arb(0.1) # TODO: How to pick this?

    (; σ) = λ

    # Precompute functions and function bounds
    p = P(ξ₁, (λ, κ))
    p_dξ = P_dξ(ξ₁, (λ, κ))
    p_dξ_dξ = P_dξ_dξ(ξ₁, (λ, κ))
    e = E(ξ₁, (λ, κ))
    e_dξ = E_dξ(ξ₁, (λ, κ))
    j_e = J_E(ξ₁, (λ, κ))
    j_p = J_P(ξ₁, (λ, κ))

    C = FunctionBounds(κ, ξ₁, λ)

    # Bounds norms
    norm_u = norm_bound_u(γ, κ, ξ₁, v, λ, C)
    norm_u_dξ = norm_bound_u_dξ(γ, κ, ξ₁, v, norm_u, λ, C)
    norm_u_dξ_dξ = norm_bound_u_dξ_dξ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ, C)
    norm_u_dξ_dξ_dξ =
        norm_bound_u_dξ_dξ_dξ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ, C)

    # Compute zeroth order bounds
    Q = add_error(zero(γ), norm_u * ξ₁^(-1 / σ + v))
    dQ = add_error(zero(γ), norm_u_dξ * ξ₁^(-1 / σ + v))

    # Improve the bounds iteratively

    # IMPROVE: In practice three iterations seems to be enough to
    # saturate the convergence. But it might be better to choose this
    # dynamically.
    for _ = 1:3
        I_P = I_P_enclose(
            γ,
            κ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dξ_dξ,
            norm_u_dξ_dξ_dξ,
            Q,
            dQ,
            λ,
            C,
            p,
            p_dξ,
            p_dξ_dξ,
        )

        Q = γ * p + e * I_P

        I_E_dξ = j_e * abs(Q)^2σ * Q
        I_P_dξ = -j_p * abs(Q)^2σ * Q

        dQ = γ * p_dξ + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ
    end

    return SVector(Q, dQ)
end

function solution_infinity(γ::ComplexF64, κ::Float64, ξ₁::Float64, λ::CGLParams{Float64})
    (; σ, d) = λ

    _, _, c = _abc(κ, λ)

    p = P(ξ₁, (λ, κ))
    p_dξ = P_dξ(ξ₁, (λ, κ))
    e = E(ξ₁, (λ, κ))
    e_dξ = E_dξ(ξ₁, (λ, κ))

    # Compute first order approximation of Q
    Q_1 = γ * p

    # Compute an improved approximation of Q and dQ
    I_E = zero(γ)
    I_P = B_W(κ, λ) * exp(-c * ξ₁^2) * p * ξ₁^(d - 2) * abs(Q_1)^2σ * Q_1 / 2c

    Q = γ * p + p * I_E + e * I_P

    I_E_dξ = J_E(ξ₁, (λ, κ)) * abs(Q)^2σ * Q
    I_P_dξ = -J_P(ξ₁, (λ, κ)) * abs(Q)^2σ * Q

    dQ = γ * p_dξ + p_dξ * I_E + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ

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

    # Precompute functions and function bounds
    p = P(ξ₁, (λ, κ))
    p_dξ = P_dξ(ξ₁, (λ, κ))
    p_dξ_dξ = P_dξ_dξ(ξ₁, (λ, κ))
    p_dκ = P_dκ(ξ₁, (λ, κ))
    p_dξ_dκ = P_dξ_dκ(ξ₁, (λ, κ))
    e = E(ξ₁, (λ, κ))
    e_dξ = E_dξ(ξ₁, (λ, κ))
    e_dκ = E_dκ(ξ₁, (λ, κ))
    e_dξ_dκ = E_dξ_dκ(ξ₁, (λ, κ))
    j_e = J_E(ξ₁, (λ, κ))
    j_p = J_P(ξ₁, (λ, κ))
    j_e_dκ = J_E_dκ(ξ₁, (λ, κ))
    j_p_dκ = J_P_dκ(ξ₁, (λ, κ))
    D_ξ₁ = D(ξ₁, (λ, κ))
    D_dξ_ξ₁ = D_dξ(ξ₁, (λ, κ))

    C = FunctionBounds(κ, ξ₁, λ, include_dκ = true)

    # Bounds norms
    norm_u = norm_bound_u(γ, κ, ξ₁, v, λ, C)
    norm_u_dξ = norm_bound_u_dξ(γ, κ, ξ₁, v, norm_u, λ, C)
    norm_u_dξ_dξ = norm_bound_u_dξ_dξ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ, C)
    norm_u_dξ_dξ_dξ =
        norm_bound_u_dξ_dξ_dξ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ, C)

    norm_u_dγ = norm_bound_u_dγ(γ, κ, ξ₁, v, norm_u, λ, C)
    norm_u_dξ_dγ = norm_bound_u_dξ_dγ(γ, κ, ξ₁, v, norm_u, norm_u_dγ, λ, C)
    norm_u_dκ = norm_bound_u_dκ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ, C)
    norm_u_dξ_dκ =
        norm_bound_u_dξ_dκ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dκ, λ, C)

    # Compute zeroth order bounds
    Q = add_error(zero(γ), norm_u * ξ₁^(-1 / σ + v))
    dQ = add_error(zero(γ), norm_u_dξ * ξ₁^(-1 / σ + v))
    Q_dγ = add_error(zero(γ), norm_u_dγ * ξ₁^(-1 / σ + v))
    dQ_dγ = add_error(zero(γ), norm_u_dξ_dγ * ξ₁^(-1 / σ + v))
    Q_dκ = add_error(zero(γ), norm_u_dκ * ξ₁^(-1 / σ + v))
    dQ_dκ = add_error(zero(γ), norm_u_dξ_dκ * ξ₁^(-1 / σ + v))

    # Improve the bounds iteratively

    # IMPROVE: In practice three iterations seems to be enough to
    # saturate the convergence. But it might be better to choose this
    # dynamically.
    for _ = 1:3
        I_P = I_P_enclose(
            γ,
            κ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dξ_dξ,
            norm_u_dξ_dξ_dξ,
            Q,
            dQ,
            λ,
            C,
            p,
            p_dξ,
            p_dξ_dξ,
        )

        I_P_dγ = I_P_dγ_enclose(
            γ,
            κ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dγ,
            norm_u_dξ_dγ,
            Q,
            Q_dγ,
            λ,
            C,
            p,
        )

        I_P_dκ = I_P_dκ_enclose(
            γ,
            κ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dξ_dξ,
            norm_u_dκ,
            norm_u_dξ_dκ,
            Q,
            dQ,
            Q_dκ,
            λ,
            C,
            p,
            D_ξ₁,
            D_dξ_ξ₁,
        )

        Q = γ * p + e * I_P
        Q_dγ = p + e * I_P_dγ
        Q_dκ = γ * p_dκ + e * I_P_dκ + e_dκ * I_P

        I_E_dξ = j_e * abs(Q)^2σ * Q
        I_P_dξ = -j_p * abs(Q)^2σ * Q

        I_E_dξ_dγ =
            j_e * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dγ) * Q + abs(Q)^2 * Q_dγ)
        I_P_dξ_dγ =
            -j_p * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dγ) * Q + abs(Q)^2 * Q_dγ)

        I_E_dξ_dκ =
            j_e_dκ * abs(Q)^2σ * Q +
            j_e * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dκ) * Q + abs(Q)^2 * Q_dκ)
        I_P_dξ_dκ =
            -j_p_dκ * abs(Q)^2σ * Q -
            j_p * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dκ) * Q + abs(Q)^2 * Q_dκ)


        dQ = γ * p_dξ + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ
        dQ_dγ = p_dξ + p * I_E_dξ_dγ + e * I_P_dξ_dγ + e_dξ * I_P_dγ
        dQ_dκ =
            γ * p_dξ_dκ +
            p * I_E_dξ_dκ +
            p_dκ * I_E_dξ +
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

"""
    solution_infinity_jacobian_epsilon(γ, κ, ξ₁, λ::CGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes Jacobian w.r.t. the parameters `γ` and `ϵ` of
`[Q(ξ₁), d(Q)(ξ₁)]. The Jacobian is given by
```
[
d(Q(ξ₁), μ) d(Q(ξ₁), ϵ)
d(d(Q)(ξ₁), μ) d((Q)(ξ₁), ϵ)
]
```
where we use `d(Q, μ)` to denote the derivative of `Q` w.r.t. `μ`.
"""
function solution_infinity_jacobian_epsilon(γ::Acb, κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    v = Arb(0.1) # TODO: How to pick this?

    (; σ) = λ

    # Precompute functions and function bounds
    p = P(ξ₁, (λ, κ))
    p_dξ = P_dξ(ξ₁, (λ, κ))
    p_dξ_dξ = P_dξ_dξ(ξ₁, (λ, κ))
    p_dϵ = P_dϵ(ξ₁, (λ, κ))
    p_dξ_dϵ = P_dξ_dϵ(ξ₁, (λ, κ))
    e = E(ξ₁, (λ, κ))
    e_dξ = E_dξ(ξ₁, (λ, κ))
    e_dϵ = E_dϵ(ξ₁, (λ, κ))
    e_dξ_dϵ = E_dξ_dϵ(ξ₁, (λ, κ))
    j_e = J_E(ξ₁, (λ, κ))
    j_p = J_P(ξ₁, (λ, κ))
    j_e_dϵ = J_E_dϵ(ξ₁, (λ, κ))
    j_p_dϵ = J_P_dϵ(ξ₁, (λ, κ))
    H_ξ₁ = H(ξ₁, (λ, κ))
    H_dξ_ξ₁ = H_dξ(ξ₁, (λ, κ))

    C = FunctionBounds(κ, ξ₁, λ, include_dϵ = true)

    # Bounds norms
    norm_u = norm_bound_u(γ, κ, ξ₁, v, λ, C)
    norm_u_dξ = norm_bound_u_dξ(γ, κ, ξ₁, v, norm_u, λ, C)
    norm_u_dξ_dξ = norm_bound_u_dξ_dξ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ, C)
    norm_u_dξ_dξ_dξ =
        norm_bound_u_dξ_dξ_dξ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ, C)

    norm_u_dγ = norm_bound_u_dγ(γ, κ, ξ₁, v, norm_u, λ, C)
    norm_u_dξ_dγ = norm_bound_u_dξ_dγ(γ, κ, ξ₁, v, norm_u, norm_u_dγ, λ, C)
    norm_u_dϵ = norm_bound_u_dϵ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ, C)
    norm_u_dξ_dϵ =
        norm_bound_u_dξ_dϵ(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dϵ, λ, C)

    # Compute zeroth order bounds
    Q = add_error(zero(γ), norm_u * ξ₁^(-1 / σ + v))
    dQ = add_error(zero(γ), norm_u_dξ * ξ₁^(-1 / σ + v))
    Q_dγ = add_error(zero(γ), norm_u_dγ * ξ₁^(-1 / σ + v))
    dQ_dγ = add_error(zero(γ), norm_u_dξ_dγ * ξ₁^(-1 / σ + v))
    Q_dϵ = add_error(zero(γ), norm_u_dϵ * ξ₁^(-1 / σ + v))
    dQ_dϵ = add_error(zero(γ), norm_u_dξ_dϵ * ξ₁^(-1 / σ + v))

    # Improve the bounds iteratively

    # IMPROVE: In practice three iterations seems to be enough to
    # saturate the convergence. But it might be better to choose this
    # dynamically.
    for _ = 1:3
        I_P = I_P_enclose(
            γ,
            κ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dξ_dξ,
            norm_u_dξ_dξ_dξ,
            Q,
            dQ,
            λ,
            C,
            p,
            p_dξ,
            p_dξ_dξ,
        )

        I_P_dγ = I_P_dγ_enclose(
            γ,
            κ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dγ,
            norm_u_dξ_dγ,
            Q,
            Q_dγ,
            λ,
            C,
            p,
        )

        I_P_dϵ = I_P_dϵ_enclose(
            γ,
            κ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dξ_dξ,
            norm_u_dϵ,
            norm_u_dξ_dϵ,
            Q,
            dQ,
            Q_dϵ,
            λ,
            C,
            p,
            H_ξ₁,
            H_dξ_ξ₁,
        )

        Q = γ * p + e * I_P
        Q_dγ = p + e * I_P_dγ
        Q_dϵ = γ * p_dϵ + e * I_P_dϵ + e_dϵ * I_P

        I_E_dξ = j_e * abs(Q)^2σ * Q
        I_P_dξ = -j_p * abs(Q)^2σ * Q

        I_E_dξ_dγ =
            j_e * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dγ) * Q + abs(Q)^2 * Q_dγ)
        I_P_dξ_dγ =
            -j_p * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dγ) * Q + abs(Q)^2 * Q_dγ)

        I_E_dξ_dϵ =
            j_e_dϵ * abs(Q)^2σ * Q +
            j_e * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dϵ) * Q + abs(Q)^2 * Q_dϵ)
        I_P_dξ_dϵ =
            -j_p_dϵ * abs(Q)^2σ * Q -
            j_p * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dϵ) * Q + abs(Q)^2 * Q_dϵ)


        dQ = γ * p_dξ + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ
        dQ_dγ = p_dξ + p * I_E_dξ_dγ + e * I_P_dξ_dγ + e_dξ * I_P_dγ
        dQ_dϵ =
            γ * p_dξ_dϵ +
            p * I_E_dξ_dϵ +
            p_dϵ * I_E_dξ +
            e * I_P_dξ_dϵ +
            e_dξ * I_P_dϵ +
            e_dϵ * I_P_dξ +
            e_dξ_dϵ * I_P
    end

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dϵ, dQ_dϵ)
end

function solution_infinity_jacobian_epsilon(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    # IMPROVE: Add higher order versions
    Q_dγ = P(ξ₁, (λ, κ))
    dQ_dγ = P_dξ(ξ₁, (λ, κ))

    Q_dϵ = γ * P_dϵ(ξ₁, (λ, κ))
    dQ_dϵ = γ * P_dξ_dϵ(ξ₁, (λ, κ))

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dϵ, dQ_dϵ)
end
