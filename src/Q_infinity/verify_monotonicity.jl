"""
    verify_monotonicity_infinity(γ, κ, ϵ, ξ₁, Λ; return_coefficients, verbose)

Return true if `abs2(Q)` could be verified to be monotone on ``[ξ₁,
∞)`` and false otherwise.

This is based on the lemmas in Section REF(sec:verify-monotonicity),
specifically Lemmas REF(lemma:Q-expansion),
REF(lemma:abs2-Q-derivative), REF(lemma:P-expansion) and
REF(lemma:monotonicity-infinity).

If `return_coefficients = Val(true)`, then also return the values of

```
lhs = abs(real(2a * abs(c^-a)^2 * (S_1 * S_2)(Arb((0, ξ₁^-1))) + R_P_conj_P_dξ))
rhs = C_R_mon / abs(p_Q)^2 * ξ₁^((2σ + 1) * v - 2)
```

corresponding to the left and right-hand sides of the condition in
REF(lemma:monotonicity-infinity).
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

    # Precompute function and norm bounds
    CU = UBounds(_abc(κ, ϵ, Λ)..., ξ₁)
    C = FunctionBounds(κ, ϵ, ξ₁, Λ, CU)
    CI = IBounds(κ, ϵ, ξ₁, v, Λ, C)

    norms = NormBounds(γ, κ, ϵ, ξ₁, v, Λ, C, CI)

    # Compute C_p_Q, p_Q, C_R_Q and C_R_dQ from Lemma REF(lemma:Q-expansion)
    C_p_Q = CI.I_E * norms.Q^(2σ + 1) * ξ₁^((2σ + 1) * v - 2)
    p_Q = add_error(γ, C_p_Q)# Enclosure of p_Q

    C_R_Q =
        C.E *
        (
            CI.I_P_2_1 * norms.Q^2 +
            CI.I_P_2_2 * norms.Q * norms.Q_dξ * ξ₁^-1 +
            CI.I_P_2_3 * norms.Q_dξ^2 +
            CI.I_P_2_4 * norms.Q * norms.Q_dξ_dξ
        ) *
        norms.Q^(2σ - 1)

    C_R_dQ =
        C.E_dξ *
        (
            CI.I_P_2_1 * norms.Q^2 +
            CI.I_P_2_2 * norms.Q * norms.Q_dξ * ξ₁^-1 +
            CI.I_P_2_3 * norms.Q_dξ^2 +
            CI.I_P_2_4 * norms.Q * norms.Q_dξ_dξ
        ) *
        norms.Q^(2σ - 1)

    # Compute C_R_mon from Lemma REF(lemma:abs2-Q-derivative)
    C_R_mon =
        (abs(γ) + C_p_Q) * C.P * C_R_dQ +
        (abs(γ) + C_p_Q) * C.P_dξ * C_R_Q * ξ₁^-2 +
        C_R_Q * C_R_dQ * ξ₁^((2σ + 1) * v - 4)

    # Compute S_1, S_2, C_R_P_conj_P_dξ and R_P_conj_P_dξ from Lemma REF(lemma:P-expansion)
    n = 10 # Number of terms in expansion

    # Compute S_1 and S_2 as polynomials in ξ^-1
    S_1 = AcbPoly()
    S_2 = AcbPoly()
    for k = 0:(n-1)
        S_1[2k] = rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k)
        S_2[2k] = conj(rising(a + 1, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k))
    end

    # Compute bounds for the sums in C_R_P_conj_P_dξ
    S_1_bound = sum(
        k ->
            abs(rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k)) * ξ₁^(-2k),
        0:(n-1),
    )
    S_2_bound = sum(
        k ->
            abs(rising(a + 1, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k)) *
            ξ₁^(-2k),
        0:(n-1),
    )

    C_R_P_conj_P_dξ =
        S_1_bound * C_R_U(n, a + 1, b + 1, c * ξ₁^2) * abs(c^-n) * ξ₁^(-2n) +
        S_2_bound * C_R_U(n, a, b, c * ξ₁^2) * abs(c^-n) * ξ₁^(-2n) +
        C_R_U(n, a, b, c * ξ₁^2) * C_R_U(n, a + 1, b + 1, c * ξ₁^2) * abs(c^-n)^2 * ξ₁^(-4n)

    # Compute an enclosure of R_P_conj_P_dξ
    R_P_conj_P_dξ = add_error(Acb(0), C_R_P_conj_P_dξ)

    # Verify the condition of Lemma REF(lemma:monotonicity-infinity)

    # Compute an enclosure of the left-hand side
    # Note that we multiple S_1 and S_2 and then evaluate them on the interval [0, ξ₁^-1]
    lhs = abs(real(2conj(a) * abs(c^-a)^2 * (S_1 * S_2)(Arb((0, ξ₁^-1))) + R_P_conj_P_dξ))
    # Compute an enclosure of the right-hand side
    rhs = C_R_mon / abs(p_Q)^2 * ξ₁^((2σ + 1) * v - 2)
    condition = lhs > rhs

    if verbose
        @info "Computed values" getinterval(lhs) rhs
    end

    if return_coefficients isa Val{false}
        return condition
    else
        return condition, lhs, rhs
    end
end
