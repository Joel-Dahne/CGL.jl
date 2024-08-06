"""
verify_monotonicity_infinity(γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}; verbose)

Verify monotonicity for `ξ >= ξ₁`.
"""
function verify_monotonicity_infinity(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    verbose = false,
)
    v = Arb(0.1) # TODO: How to pick this?

    (; σ) = λ

    a, b, c = _abc(κ, ϵ, λ)

    # Precompute functions as well as function and norm bounds
    F = FunctionEnclosures(ξ₁, κ, ϵ, λ)

    CU = UBounds(_abc(κ, ϵ, λ)..., ξ₁)
    C = FunctionBounds(κ, ϵ, ξ₁, λ, CU)

    norms = NormBounds(γ, κ, ϵ, ξ₁, v, λ, C)

    # Compute needed bounds
    p_Q = add_error(
        c^-a * γ,
        abs(c^-a) * C_I_E(κ, ϵ, ξ₁, v, λ, C) * norms.Q^(2σ + 1) * ξ₁^((2σ + 1) * v - 2),
    )

    abs_p_X = 4abs2(p_Q) * abs(real(a))

    R_Q_bound =
        C_R_U(1, a, b, c * ξ₁^2) * abs(c^-1 * p_Q) * ξ₁^((-2σ + 1) * v) +
        C.E * C_I_P(κ, ϵ, ξ₁, v, λ, C) * norms.Q^(2σ + 1)

    R_dQ_bound =
        C_R_U(1, a + 1, b + 1, c * ξ₁^2) * abs(2a * c^-1 * p_Q) * ξ₁^((-2σ + 1) * v) +
        C.P * C.J_E * norms.Q^(2σ + 1) +
        C.E_dξ *
        (
            C_I_P_2_1(κ, ϵ, ξ₁, v, λ, C) * norms.Q^2 +
            C_I_P_2_2(κ, ϵ, ξ₁, v, λ, C) * norms.Q^2 * ξ₁^-2 +
            C_I_P_2_3(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ * ξ₁^-1 +
            C_I_P_2_4(κ, ϵ, ξ₁, v, λ, C) * norms.Q_dξ^2 +
            C_I_P_2_5(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ_dξ
        ) *
        norms.Q^(2σ - 1) +
        C.E * C.J_P * norms.Q^(2σ + 1)

    R_X_bound =
        4abs(p_Q * R_dQ_bound) +
        8abs(a * p_Q * R_Q_bound) * ξ₁^-1 +
        4abs(R_Q_bound * R_dQ_bound) * ξ₁^((2σ + 1) * v - 3)

    if verbose
        @info "Computed values" abs_p_X R_X_bound R_X_bound * ξ₁^((2σ + 1) * v - 2)
    end

    return R_X_bound * ξ₁^((2σ + 1) * v - 2) <= abs_p_X
end
