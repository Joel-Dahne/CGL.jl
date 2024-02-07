"""
    solution_infinity(γ, κ, ϵ, ξ₁, λ::CGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes `[Q(ξ₁), d(Q)(ξ₁)]`.
"""
function solution_infinity(γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    v = Arb(0.1) # TODO: How to pick this?

    (; σ) = λ

    # Precompute functions and function bounds
    F = FunctionEnclosures(ξ₁, κ, ϵ, λ)
    C = FunctionBounds(κ, ϵ, ξ₁, λ)

    # Bounds norms
    norm_u = norm_bound_u(γ, κ, ϵ, ξ₁, v, λ, C)
    norm_u_dξ = norm_bound_u_dξ(γ, κ, ϵ, ξ₁, v, norm_u, λ, C)
    norm_u_dξ_dξ = norm_bound_u_dξ_dξ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dξ, λ, C)
    norm_u_dξ_dξ_dξ =
        norm_bound_u_dξ_dξ_dξ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ, C)

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
            ϵ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dξ_dξ,
            norm_u_dξ_dξ_dξ,
            Q,
            dQ,
            λ,
            F,
            C,
        )

        Q = γ * F.P + F.E * I_P

        I_E_dξ = F.J_E * abs(Q)^2σ * Q
        I_P_dξ = -F.J_P * abs(Q)^2σ * Q

        dQ = γ * F.P_dξ + F.P * I_E_dξ + F.E_dξ * I_P + F.E * I_P_dξ
    end

    return SVector(Q, dQ)
end

function solution_infinity(
    γ::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, ϵ, λ)

    # Precompute functions
    F = FunctionEnclosures(ξ₁, κ, ϵ, λ)

    # Compute first order approximation of Q
    Q_1 = γ * F.P

    # Compute an improved approximation of Q and dQ
    I_E = zero(γ)
    I_P = B_W(κ, ϵ, λ) * exp(-c * ξ₁^2) * F.P * ξ₁^(d - 2) * abs(Q_1)^2σ * Q_1 / 2c

    Q = γ * F.P + F.P * I_E + F.E * I_P

    I_E_dξ = F.J_E * abs(Q)^2σ * Q
    I_P_dξ = -F.J_P * abs(Q)^2σ * Q

    dQ = γ * F.P_dξ + F.P_dξ * I_E + F.P * I_E_dξ + F.E_dξ * I_P + F.E * I_P_dξ

    return SVector(Q, dQ)
end

"""
    solution_infinity_jacobian_kappa(γ, κ, ϵ, ξ₁, λ::CGLParams)

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
function solution_infinity_jacobian_kappa(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    v = Arb(0.1) # TODO: How to pick this?

    (; σ) = λ

    # Precompute functions and function bounds
    F = FunctionEnclosures(ξ₁, κ, ϵ, λ, include_dκ = true)
    C = FunctionBounds(κ, ϵ, ξ₁, λ, include_dκ = true)

    # Bounds norms
    norm_u = norm_bound_u(γ, κ, ϵ, ξ₁, v, λ, C)
    norm_u_dξ = norm_bound_u_dξ(γ, κ, ϵ, ξ₁, v, norm_u, λ, C)
    norm_u_dξ_dξ = norm_bound_u_dξ_dξ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dξ, λ, C)
    norm_u_dξ_dξ_dξ =
        norm_bound_u_dξ_dξ_dξ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ, C)

    norm_u_dγ = norm_bound_u_dγ(γ, κ, ϵ, ξ₁, v, norm_u, λ, C)
    norm_u_dξ_dγ = norm_bound_u_dξ_dγ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dγ, λ, C)
    norm_u_dκ = norm_bound_u_dκ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ, C)
    norm_u_dξ_dκ =
        norm_bound_u_dξ_dκ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dκ, λ, C)

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
            ϵ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dξ_dξ,
            norm_u_dξ_dξ_dξ,
            Q,
            dQ,
            λ,
            F,
            C,
        )

        I_P_dγ = I_P_dγ_enclose(
            γ,
            κ,
            ϵ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dγ,
            norm_u_dξ_dγ,
            Q,
            Q_dγ,
            λ,
            F,
            C,
        )

        I_P_dκ = I_P_dκ_enclose(
            γ,
            κ,
            ϵ,
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
            F,
            C,
        )

        Q = γ * F.P + F.E * I_P
        Q_dγ = F.P + F.E * I_P_dγ
        Q_dκ = γ * F.P_dκ + F.E * I_P_dκ + F.E_dκ * I_P

        I_E_dξ = F.J_E * abs(Q)^2σ * Q
        I_P_dξ = -F.J_P * abs(Q)^2σ * Q

        I_E_dξ_dγ =
            F.J_E * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dγ) * Q + abs(Q)^2 * Q_dγ)
        I_P_dξ_dγ =
            -F.J_P * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dγ) * Q + abs(Q)^2 * Q_dγ)

        I_E_dξ_dκ =
            F.J_E_dκ * abs(Q)^2σ * Q +
            F.J_E * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dκ) * Q + abs(Q)^2 * Q_dκ)
        I_P_dξ_dκ =
            -F.J_P_dκ * abs(Q)^2σ * Q -
            F.J_P * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dκ) * Q + abs(Q)^2 * Q_dκ)


        dQ = γ * F.P_dξ + F.P * I_E_dξ + F.E_dξ * I_P + F.E * I_P_dξ
        dQ_dγ = F.P_dξ + F.P * I_E_dξ_dγ + F.E * I_P_dξ_dγ + F.E_dξ * I_P_dγ
        dQ_dκ =
            γ * F.P_dξ_dκ +
            F.P * I_E_dξ_dκ +
            F.P_dκ * I_E_dξ +
            F.E * I_P_dξ_dκ +
            F.E_dξ * I_P_dκ +
            F.E_dκ * I_P_dξ +
            F.E_dξ_dκ * I_P
    end

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end

function solution_infinity_jacobian_kappa(
    γ::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    # Precompute functions
    F = FunctionEnclosures(ξ₁, κ, ϵ, λ, include_dκ = true)

    # IMPROVE: Add higher order versions
    Q_dγ = F.P
    dQ_dγ = F.P_dξ

    Q_dκ = γ * F.P_dκ
    dQ_dκ = γ * F.P_dξ_dκ

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end

"""
    solution_infinity_jacobian_epsilon(γ, κ, ϵ, ξ₁, λ::CGLParams)

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
function solution_infinity_jacobian_epsilon(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    v = Arb(0.1) # TODO: How to pick this?

    (; σ) = λ

    # Precompute functions and function bounds
    F = FunctionEnclosures(ξ₁, κ, ϵ, λ, include_dϵ = true)
    C = FunctionBounds(κ, ϵ, ξ₁, λ, include_dϵ = true)

    # Bounds norms
    norm_u = norm_bound_u(γ, κ, ϵ, ξ₁, v, λ, C)
    norm_u_dξ = norm_bound_u_dξ(γ, κ, ϵ, ξ₁, v, norm_u, λ, C)
    norm_u_dξ_dξ = norm_bound_u_dξ_dξ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dξ, λ, C)
    norm_u_dξ_dξ_dξ =
        norm_bound_u_dξ_dξ_dξ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ, C)

    norm_u_dγ = norm_bound_u_dγ(γ, κ, ϵ, ξ₁, v, norm_u, λ, C)
    norm_u_dξ_dγ = norm_bound_u_dξ_dγ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dγ, λ, C)
    norm_u_dϵ = norm_bound_u_dϵ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ, C)
    norm_u_dξ_dϵ =
        norm_bound_u_dξ_dϵ(γ, κ, ϵ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dϵ, λ, C)

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
            ϵ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dξ_dξ,
            norm_u_dξ_dξ_dξ,
            Q,
            dQ,
            λ,
            F,
            C,
        )

        I_P_dγ = I_P_dγ_enclose(
            γ,
            κ,
            ϵ,
            ξ₁,
            v,
            norm_u,
            norm_u_dξ,
            norm_u_dγ,
            norm_u_dξ_dγ,
            Q,
            Q_dγ,
            λ,
            F,
            C,
        )

        I_P_dϵ = I_P_dϵ_enclose(
            γ,
            κ,
            ϵ,
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
            F,
            C,
        )

        Q = γ * F.P + F.E * I_P
        Q_dγ = F.P + F.E * I_P_dγ
        Q_dϵ = γ * F.P_dϵ + F.E * I_P_dϵ + F.E_dϵ * I_P

        I_E_dξ = F.J_E * abs(Q)^2σ * Q
        I_P_dξ = -F.J_P * abs(Q)^2σ * Q

        I_E_dξ_dγ =
            F.J_E * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dγ) * Q + abs(Q)^2 * Q_dγ)
        I_P_dξ_dγ =
            -F.J_P * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dγ) * Q + abs(Q)^2 * Q_dγ)

        I_E_dξ_dϵ =
            F.J_E_dϵ * abs(Q)^2σ * Q +
            F.J_E * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dϵ) * Q + abs(Q)^2 * Q_dϵ)
        I_P_dξ_dϵ =
            -F.J_P_dϵ * abs(Q)^2σ * Q -
            F.J_P * abs(Q)^(2σ - 2) * (2σ * real(conj(Q) * Q_dϵ) * Q + abs(Q)^2 * Q_dϵ)


        dQ = γ * F.P_dξ + F.P * I_E_dξ + F.E_dξ * I_P + F.E * I_P_dξ
        dQ_dγ = F.P_dξ + F.P * I_E_dξ_dγ + F.E * I_P_dξ_dγ + F.E_dξ * I_P_dγ
        dQ_dϵ =
            γ * F.P_dξ_dϵ +
            F.P * I_E_dξ_dϵ +
            F.P_dϵ * I_E_dξ +
            F.E * I_P_dξ_dϵ +
            F.E_dξ * I_P_dϵ +
            F.E_dϵ * I_P_dξ +
            F.E_dξ_dϵ * I_P
    end

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dϵ, dQ_dϵ)
end

function solution_infinity_jacobian_epsilon(
    γ::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    # Precompute functions
    F = FunctionEnclosures(ξ₁, κ, ϵ, λ, include_dϵ = true)

    # IMPROVE: Add higher order versions
    Q_dγ = F.P
    dQ_dγ = F.P_dξ

    Q_dϵ = γ * F.P_dϵ
    dQ_dϵ = γ * F.P_dξ_dϵ

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dϵ, dQ_dϵ)
end
