"""
    verify_monotonicity_infinity(γ, κ, ϵ, ξ₁, Λ; return_coefficients, verbose)

Return `ξ₂` such that `abs2(Q)` is monotone on ``[ξ₂, ∞)``. If no such
`ξ₂` could be found, return an indeterminate value.

This is based on Lemmas REF(lemma:Q-expansion) and
REF(lemma:monotonicity-infinity).

If `return_coefficients = Val(true)`, then also return the values of
`C_p_mon, C_R_mon, ξ₂_lower_bound` that are used in the computation of
`ξ₂`.
"""
function verify_monotonicity_infinity(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb};
    return_coefficients::Union{Val{false},Val{true}} = Val{false}(),
    verbose = false,
)
    v = Arb("0.1")

    (; σ) = Λ

    a, b, c = _abc(κ, ϵ, Λ)

    # Precompute functions as well as function and norm bounds
    CU = UBounds(_abc(κ, ϵ, Λ)..., ξ₁)
    C = FunctionBounds(κ, ϵ, ξ₁, Λ, CU)
    CI = IBounds(κ, ϵ, ξ₁, v, Λ, C)

    norms = NormBounds(γ, κ, ϵ, ξ₁, v, Λ, C, CI)

    # Compute C_p_Q, C_R_Q and C_R_dQ from Lemma REF(lemma:Q-expansion)
    n = 10 # Number of terms in expansion when bounding P and P_dξ

    C_p_Q = abs(c^-a) * CI.I_E * norms.Q^(2σ + 1) * ξ₁^((2σ + 1) * v - 2)

    C_R_Q =
        (abs(c^-a * γ) + C_p_Q) *
        (
            sum(
                k ->
                    abs(rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k)) *
                    ξ₁^(-2k + 2),
                1:(n-1),
            ) + C_R_U(n, a, b, c * ξ₁^2) * abs(c^-n) * ξ₁^(-2n + 2)
        ) *
        ξ₁^((-2σ + 1) * v) + C.E * CI.I_P * norms.Q^(2σ + 1)

    C_R_dQ =
        abs(2a) *
        (abs(c^-a * γ) + C_p_Q) *
        (
            sum(
                k ->
                    abs(rising(a + 1, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k)) *
                    ξ₁^(-2k + 2),
                1:(n-1),
            ) + C_R_U(n, a + 1, b + 1, c * ξ₁^2) * abs(c^-n) * ξ₁^(-2n + 2)
        ) *
        ξ₁^((-2σ + 1) * v) +
        C.P * C.J_E * norms.Q^(2σ + 1) +
        C.E_dξ *
        (
            CI.I_P_2_1 * norms.Q^2 +
            CI.I_P_2_3 * norms.Q * norms.Q_dξ * ξ₁^-1 +
            CI.I_P_2_4 * norms.Q_dξ^2 +
            CI.I_P_2_5 * norms.Q * norms.Q_dξ_dξ
        ) *
        norms.Q^(2σ - 1) +
        C.E * C.J_P * norms.Q^(2σ + 1)

    # Compute C_p_mon and C_R_mon from Lemma REF(lemma:monotonicity-infinity)
    C_p_mon = if abs(c^-a * γ) > C_p_Q
        4abs(real(a)) * (abs(c^-a * γ) - C_p_Q)^2
    else
        indeterminate(Arb)
    end

    C_R_mon =
        4C_p_Q * C_R_dQ +
        8abs(a) * C_p_Q * C_R_Q * ξ₁^-1 +
        4C_R_Q * C_R_dQ * ξ₁^((2σ + 1) * v - 3)

    # For ξ greater than this we have that abs2(Q) is monotone.
    ξ₂_lower_bound = (C_p_mon / C_R_mon)^inv((2σ + 1) * v - 2)
    # We can take a lower bound than ξ₁
    ξ₂ = max(ξ₁, ξ₂_lower_bound)

    if verbose
        @info "Computed values" C_p_Q C_R_Q C_R_dQ C_p_mon C_R_mon ξ₂_lower_bound
    end

    if return_coefficients isa Val{false}
        return ξ₂
    else
        return ξ₂, C_p_mon, C_R_mon, ξ₂_lower_bound
    end
end
