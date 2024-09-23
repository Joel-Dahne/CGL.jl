"""
verify_monotonicity_infinity(γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}; verbose)

Return `ξ₂` such that `abs2(Q)` is monotone on ``[ξ₂, ∞)``. If no such
`ξ₂` could be found, return and indeterminate value.
"""
function verify_monotonicity_infinity(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    return_coefficients::Union{Val{false},Val{true}} = Val{false}(),
    verbose = false,
)
    v = Arb(0.1) # TODO: How to pick this?

    (; σ) = λ

    a, b, c = _abc(κ, ϵ, λ)

    # Precompute functions as well as function and norm bounds
    CU = UBounds(_abc(κ, ϵ, λ)..., ξ₁)
    C = FunctionBounds(κ, ϵ, ξ₁, λ, CU)

    norms = NormBounds(γ, κ, ϵ, ξ₁, v, λ, C)

    # Compute needed bounds
    n = 10 # Number of terms in expansion when bounding P and P_dξ

    C_p_Q = abs(c^-a) * C_I_E(κ, ϵ, ξ₁, v, λ, C) * norms.Q^(2σ + 1) * ξ₁^((2σ + 1) * v - 2)

    C_R_Q =
        (abs(c^-a * γ) + C_p_Q) *
        (
            sum(
                k ->
                    abs(rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k)) *
                    ξ₁^(-2k + 2),
                1:n-1,
            ) + C_R_U(n, a, b, c * ξ₁^2) * abs(c^-n) * ξ₁^(-2n + 2)
        ) *
        ξ₁^((-2σ + 1) * v) + C.E * C_I_P(κ, ϵ, ξ₁, v, λ, C) * norms.Q^(2σ + 1)

    C_R_dQ =
        abs(2a) *
        (abs(c^-a * γ) + C_p_Q) *
        (
            sum(
                k ->
                    abs(rising(a + 1, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k)) *
                    ξ₁^(-2k + 2),
                1:n-1,
            ) + C_R_U(n, a + 1, b + 1, c * ξ₁^2) * abs(c^-n) * ξ₁^(-2n + 2)
        ) *
        ξ₁^((-2σ + 1) * v) +
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

    C_R_X =
        4C_p_Q * C_R_dQ +
        8abs(a) * C_p_Q * C_R_Q * ξ₁^-1 +
        4C_R_Q * C_R_dQ * ξ₁^((2σ + 1) * v - 3)

    abs_p_X_lower = 4abs(real(a)) * (abs(c^-a * γ) - C_p_Q)^2

    if verbose
        @info "Computed values" C_p_Q C_R_Q C_R_dQ C_R_X abs_p_X_lower
    end

    ξ₂ = if !Arblib.ispositive(abs_p_X_lower)
        verbose && @warn "Non-positive bound lower bound for abs(p_X)" abs_p_X_lower

        indeterminate(Arb)
    else
        max(ξ₁, (abs_p_X_lower / C_R_X)^inv((2σ + 1) * v - 2))
    end

    if return_coefficients isa Val{false}
        return ξ₂
    else
        return ξ₂, C_R_X, abs_p_X_lower
    end
end
