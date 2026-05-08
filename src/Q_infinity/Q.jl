"""
    Q_infinity(γ, κ, ϵ, ξ₁, Λ::CGLParams)

Compute the solution to the ODE on the interval ``[ξ₁, ∞)``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Q_infinity(γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb})
    v = Arb("0.1")

    (; σ) = Λ

    # Precompute functions as well as function and norm bounds
    F = FunctionEnclosures(κ, ϵ, ξ₁, Λ)

    CU = UBounds(_abc(κ, ϵ, Λ)..., ξ₁)
    C = FunctionBounds(κ, ϵ, ξ₁, Λ, CU)
    CI = IBounds(κ, ϵ, ξ₁, v, Λ, C)

    norms = NormBounds(γ, κ, ϵ, ξ₁, v, Λ, C, CI)

    # Compute zeroth order bounds
    Q = add_error(zero(γ), norms.Q * ξ₁^(-1 / σ + v))
    dQ = add_error(zero(γ), norms.Q_dξ * ξ₁^(-1 / σ + v))

    # Improve the bounds iteratively. Three iterations seems to be
    # enough to saturate it.
    for _ = 1:3
        I_P = I_P_enclose(γ, κ, ϵ, ξ₁, v, Q, dQ, Λ, F, C, norms)

        Q = γ * F.P + F.E * I_P

        I_E_dξ = F.J_E * abs(Q)^2σ * Q
        I_P_dξ = -F.J_P * abs(Q)^2σ * Q

        dQ = γ * F.P_dξ + F.P * I_E_dξ + F.E_dξ * I_P + F.E * I_P_dξ
    end

    return SVector(Q, dQ)
end

function Q_infinity(
    γ::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    Λ::CGLParams{Float64},
)
    (; d, σ) = Λ

    _, _, c = _abc(κ, ϵ, Λ)

    # Precompute functions
    p = P(ξ₁, κ, ϵ, Λ)
    p_dξ = P_dξ(ξ₁, κ, ϵ, Λ)
    e = E(ξ₁, κ, ϵ, Λ)
    e_dξ = E_dξ(ξ₁, κ, ϵ, Λ)
    j_p = J_P(ξ₁, κ, ϵ, Λ)
    j_e = J_E(ξ₁, κ, ϵ, Λ)

    # Compute first order approximation of Q
    Q = γ * p

    # Compute an improved approximation of Q and dQ
    I_P = B_W(κ, ϵ, Λ) * exp(-c * ξ₁^2) * p * ξ₁^(d - 2) * abs(Q)^2σ * Q / 2c

    Q = γ * p + e * I_P

    I_E_dξ = j_e * abs(Q)^2σ * Q
    I_P_dξ = -j_p * abs(Q)^2σ * Q

    dQ = γ * p_dξ + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ

    return SVector(Q, dQ)
end

"""
    Q_infinity_jacobian_kappa(γ, κ, ϵ, ξ₁, Λ::CGLParams)

This function computes the Jacobian of [`Q_infinity`](@ref) w.r.t. the
parameters `γ` and `κ`.
"""
function Q_infinity_jacobian_kappa(γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb})
    v = Arb("0.1")

    (; σ) = Λ

    # Precompute functions as well as function and norm bounds
    F = FunctionEnclosures(κ, ϵ, ξ₁, Λ, include_dκ = true)

    CU = UBounds(_abc(κ, ϵ, Λ)..., ξ₁, include_da = true)
    C = FunctionBounds(κ, ϵ, ξ₁, Λ, CU, include_dκ = true)
    CI = IBounds(κ, ϵ, ξ₁, v, Λ, C, include_dκ = true)

    norms = NormBounds(γ, κ, ϵ, ξ₁, v, Λ, C, CI, include_dκ = true)

    # Compute zeroth order bounds
    Q = add_error(zero(γ), norms.Q * ξ₁^(-1 / σ + v))
    dQ = add_error(zero(γ), norms.Q_dξ * ξ₁^(-1 / σ + v))
    Q_dγ = add_error(zero(γ), norms.Q_dγ * ξ₁^(-1 / σ + v))
    dQ_dγ = add_error(zero(γ), norms.Q_dγ_dξ * ξ₁^(-1 / σ + v))
    Q_dκ = add_error(zero(γ), norms.Q_dκ * ξ₁^(-1 / σ + v))
    dQ_dκ = add_error(zero(γ), norms.Q_dκ_dξ * ξ₁^(-1 / σ + v))

    # Improve the bounds iteratively. Three iterations seems to be
    # enough to saturate it.
    for _ = 1:3
        I_P = I_P_enclose(γ, κ, ϵ, ξ₁, v, Q, dQ, Λ, F, C, norms)

        I_P_dγ = I_P_dγ_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dγ, Λ, F, C, norms)

        I_P_dκ = I_P_dκ_enclose(γ, κ, ϵ, ξ₁, v, Q, dQ, Q_dκ, Λ, F, C, norms)

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

function Q_infinity_jacobian_kappa(
    γ::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    Λ::CGLParams{Float64},
)
    # IMPROVE: Add higher order versions
    Q_dγ = P(ξ₁, κ, ϵ, Λ)
    dQ_dγ = P_dξ(ξ₁, κ, ϵ, Λ)

    Q_dκ = γ * P_dκ(ξ₁, κ, ϵ, Λ)
    dQ_dκ = γ * P_dξ_dκ(ξ₁, κ, ϵ, Λ)

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end

"""
    Q_infinity_jacobian_epsilon(γ, κ, ϵ, ξ₁, Λ::CGLParams)

This function computes the Jacobian of [`Q_infinity`](@ref) w.r.t. the
parameters `μ` and `ϵ`.
"""
function Q_infinity_jacobian_epsilon(γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb})
    v = Arb("0.1")

    (; σ) = Λ

    # Precompute functions as well as function and norm bounds
    F = FunctionEnclosures(κ, ϵ, ξ₁, Λ, include_dϵ = true)

    CU = UBounds(_abc(κ, ϵ, Λ)..., ξ₁, include_da = true)
    C = FunctionBounds(κ, ϵ, ξ₁, Λ, CU, include_dϵ = true)
    CI = IBounds(κ, ϵ, ξ₁, v, Λ, C, include_dϵ = true)

    norms = NormBounds(γ, κ, ϵ, ξ₁, v, Λ, C, CI, include_dϵ = true)

    # Compute zeroth order bounds
    Q = add_error(zero(γ), norms.Q * ξ₁^(-1 / σ + v))
    dQ = add_error(zero(γ), norms.Q_dξ * ξ₁^(-1 / σ + v))
    Q_dγ = add_error(zero(γ), norms.Q_dγ * ξ₁^(-1 / σ + v))
    dQ_dγ = add_error(zero(γ), norms.Q_dγ_dξ * ξ₁^(-1 / σ + v))
    Q_dϵ = add_error(zero(γ), norms.Q_dϵ * ξ₁^(-1 / σ + v))
    dQ_dϵ = add_error(zero(γ), norms.Q_dϵ_dξ * ξ₁^(-1 / σ + v))

    # Improve the bounds iteratively. Three iterations seems to be
    # enough to saturate it.
    for _ = 1:3
        I_P = I_P_enclose(γ, κ, ϵ, ξ₁, v, Q, dQ, Λ, F, C, norms)

        I_P_dγ = I_P_dγ_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dγ, Λ, F, C, norms)

        I_P_dϵ = I_P_dϵ_enclose(γ, κ, ϵ, ξ₁, v, Q, dQ, Q_dϵ, Λ, F, C, norms)

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

function Q_infinity_jacobian_epsilon(
    γ::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    Λ::CGLParams{Float64},
)
    # IMPROVE: Add higher order versions
    Q_dγ = P(ξ₁, κ, ϵ, Λ)
    dQ_dγ = P_dξ(ξ₁, κ, ϵ, Λ)

    Q_dϵ = γ * P_dϵ(ξ₁, κ, ϵ, Λ)
    dQ_dϵ = γ * P_dξ_dϵ(ξ₁, κ, ϵ, Λ)

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dϵ, dQ_dϵ)
end
