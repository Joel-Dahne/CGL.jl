"""
    I_P_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dξ, Λ, F, C, norms)

Compute an enclosure of ``I_P`` using the expansion from Lemma
REF(lemma:I_P-expansion) with `n = 3`.
"""
function I_P_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Q::Acb,
    Q_dξ::Acb,
    Λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, ω, σ, δ) = Λ
    c = _c(κ, ϵ, Λ)

    # Compute enclosures of abs(Q)^2σ * Q and its first two
    # derivatives at ξ = ξ₁
    Q2σQ, Q2σQ_dξ, Q2σQ_dξ_dξ = let
        # Enclose Q_dξ_dξ using both the differential equation and the
        # bound for the norm and take the intersection.
        Q_dξ_dξ_1 = -(d - 1) / ξ₁ * Q_dξ
        -(im * κ * ξ₁ * Q_dξ + im * κ / σ * Q - ω * Q + (1 + im * δ) * abs(Q)^2σ * Q) /
        (1 - im * ϵ)

        Q_dξ_dξ_2 = add_error(zero(γ), norms.Q_dξ_dξ * ξ₁^(-1 / σ + v))

        Q_dξ_dξ = Acb(
            Arblib.intersection(real(Q_dξ_dξ_1), real(Q_dξ_dξ_2)),
            Arblib.intersection(imag(Q_dξ_dξ_1), imag(Q_dξ_dξ_2)),
        )

        a = ArbSeries((real(Q), real(Q_dξ), real(Q_dξ_dξ)))
        b = ArbSeries((imag(Q), imag(Q_dξ), imag(Q_dξ_dξ)))

        Q2σQ = abspow(a^2 + b^2, σ) * (a + im * b)

        Q2σQ[0], Q2σQ[1], 2Q2σQ[2]
    end

    # For computing I_P_2 and I_P_3 we explicitly expand the
    # derivatives in the formulas in Lemma REF(lemma:I_P-expansion)

    I_P_1 = exp(-c * ξ₁^2) * F.P * ξ₁^(d - 2) * Q2σQ

    I_P_2 =
        exp(-c * ξ₁^2) * (
            F.P_dξ * ξ₁^(d - 3) * Q2σQ +
            (d - 2) * F.P * ξ₁^(d - 4) * Q2σQ +
            F.P * ξ₁^(d - 3) * Q2σQ_dξ
        )

    I_P_3 =
        exp(-c * ξ₁^2) * (
            F.P_dξ_dξ * ξ₁^(d - 4) * Q2σQ +
            (2d - 5) * F.P_dξ * ξ₁^(d - 5) * Q2σQ +
            2F.P_dξ * ξ₁^(d - 4) * Q2σQ_dξ +
            (d - 2) * (d - 4) * F.P * ξ₁^(d - 6) * Q2σQ +
            (2d - 5) * F.P * ξ₁^(d - 5) * Q2σQ_dξ +
            F.P * ξ₁^(d - 5) * Q2σQ_dξ_dξ
        )

    # Compute bound of hat_I_P_4. Based on Lemma REF(lemma:I_P-remainder-bounds).
    α = (2σ + 1) * v - 2 / σ + d # The lemma uses α to denote this value
    @assert α - 2 < 0 # Requirement for lemma

    hat_I_P_4_bound =
        (
            (
                C.P_dξ_dξ_dξ +
                abs(3d - 9) * C.P_dξ_dξ +
                abs(3d^2 - 21d + 33) * C.P_dξ +
                abs((d - 2) * (d - 4) * (d - 6)) * C.P
            ) / abs(α - 8) *
            norms.Q^(2σ + 1) *
            ξ₁^-3 +
            (2σ + 1) * (3C.P_dξ_dξ + abs(6d - 18) * C.P_dξ + abs(3d^2 - 21d + 33) * C.P) /
            abs(α - 7) *
            norms.Q^2σ *
            norms.Q_dξ *
            ξ₁^-2 +
            (2σ + 1) * (3C.P_dξ + abs(3d - 9) * C.P) / abs(α - 6) *
            (2σ * norms.Q_dξ^2 + norms.Q * norms.Q_dξ_dξ) *
            norms.Q^(2σ - 1) *
            ξ₁^-1 +
            (2σ + 1) * C.P / abs(α - 5) *
            (
                2σ * (2σ - 1) * norms.Q_dξ^3 +
                6σ * norms.Q * norms.Q_dξ * norms.Q_dξ_dξ +
                norms.Q^2 * norms.Q_dξ_dξ_dξ
            ) *
            norms.Q^(2σ - 2)
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^(α - 5)

    main = B_W(κ, ϵ, Λ) * (I_P_1 / 2c + I_P_2 / (2c)^2 + I_P_3 / (2c)^3)
    remainder = add_error(zero(γ), abs(B_W(κ, ϵ, Λ) / (2c)^3) * hat_I_P_4_bound)

    return main + remainder
end

"""
    I_P_dγ_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dγ, Λ, F, C, norms)

Compute an enclosure of ``I_P_dγ`` using the expansion from Lemma
REF(lemma:I_P-derivatives-expansion-1).
"""
function I_P_dγ_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Q::Acb,
    Q_dγ::Acb,
    Λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, σ) = Λ
    c = _c(κ, ϵ, Λ)


    # Compute enclosure of abs(Q)^2σ * Q differentiated w.r.t γ at ξ = ξ₁
    Q2σQ_dγ = let
        a = ArbSeries((real(Q), real(Q_dγ)))
        b = ArbSeries((imag(Q), imag(Q_dγ)))

        Q2σQ = abspow(a^2 + b^2, σ) * (a + im * b)

        Q2σQ[1]
    end

    I_P_dγ_1 = exp(-c * ξ₁^2) * F.P * ξ₁^(d - 2) * Q2σQ_dγ

    # Compute bound of hat_I_P_dγ_2. This is based on Lemma
    # REF(lemma:I_P-remainder-bounds-1)
    α = (2σ + 1) * v - 2 / σ + d # The lemma uses α to denote this value
    @assert α - 2 < 0 # Requirement for lemma

    hat_I_P_dγ_2_bound =
        (
            (2σ + 1) * (C.P_dξ + abs(d - 2) * C.P) / abs(α - 4) *
            norms.Q^2σ *
            norms.Q_dγ *
            ξ₁^-1 +
            (2σ + 1) * C.P / abs(α - 3) *
            (2σ * norms.Q_dξ * norms.Q_dγ + norms.Q * norms.Q_dγ_dξ) *
            norms.Q^(2σ - 1)
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^(α - 3)

    main = B_W(κ, ϵ, Λ) * (I_P_dγ_1 / 2c)
    remainder = add_error(zero(γ), abs(B_W(κ, ϵ, Λ) / 2c) * hat_I_P_dγ_2_bound)

    return main + remainder
end

function I_P_dκ_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Q::Acb,
    Q_dξ::Acb,
    Q_dκ::Acb,
    Λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    return I_P_dκ_1_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dξ, Λ, F, C, norms) +
           I_P_dκ_2_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dκ, Λ, F, C, norms)
end

"""
    I_P_dκ_1_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dξ, Λ, F, C, norms)

Compute an enclosure of ``I_P_dκ_1`` using the expansion from Lemma
REF(lemma:I_P-derivatives-expansion-2).
"""
function I_P_dκ_1_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Q::Acb,
    Q_dξ::Acb,
    Λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, σ) = Λ
    c = _c(κ, ϵ, Λ)

    # Compute enclosure of abs(Q)^2σ * Q and its first derivative at ξ = ξ₁
    Q2σQ, Q2σQ_dξ = let
        a = ArbSeries((real(Q), real(Q_dξ)))
        b = ArbSeries((imag(Q), imag(Q_dξ)))

        Q2σQ = abspow(a^2 + b^2, σ) * (a + im * b)

        Q2σQ[0], Q2σQ[1]
    end

    # For computing I_P_dκ_1_2 we explicitly expand the derivative in
    # the formulas in Lemma REF(lemma:I_P-derivatives-expansion-2)

    I_P_dκ_1_1 = exp(-c * ξ₁^2) * F.D * ξ₁^d * Q2σQ

    I_P_dκ_1_2 =
        exp(-c * ξ₁^2) * (
            F.D_dξ * ξ₁^(d - 1) * Q2σQ +
            d * F.D * ξ₁^(d - 2) * Q2σQ +
            F.D * ξ₁^(d - 1) * Q2σQ_dξ
        )

    # Compute bound of hat_I_P_dκ_1_3. This is based on Lemma
    # REF(lemma:I_P-remainder-bounds-2)
    α = (2σ + 1) * v - 2 / σ + d # The lemma uses α to denote this value
    @assert α - 2 < 0 # Requirement for lemma

    hat_I_P_dκ_1_3_bound =
        (
            (C.D_dξ_dξ + abs(2d - 1) * C.D_dξ + abs(d * (d - 2)) * C.D) / abs(α - 4) *
            norms.Q^(2σ + 1) *
            ξ₁^-2 +
            (2σ + 1) * (2C.D_dξ + abs(2d - 1) * C.D) / abs(α - 3) *
            norms.Q^2σ *
            norms.Q_dξ *
            ξ₁^-1 +
            (2σ + 1) * C.D / abs(α - 2) *
            (2σ * norms.Q_dξ^2 + norms.Q * norms.Q_dξ_dξ) *
            norms.Q^(2σ - 1)
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^(α - 2)

    main = I_P_dκ_1_1 / 2c + I_P_dκ_1_2 / (2c)^2
    remainder = add_error(zero(γ), abs(1 / (2c)^2) * hat_I_P_dκ_1_3_bound)

    return main + remainder
end

"""
    I_P_dκ_2_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dκ, Λ, F, C, norms)

Compute an enclosure of ``I_P_dκ_2`` using the expansion from Lemma
REF(lemma:I_P-derivatives-expansion-1).
"""
function I_P_dκ_2_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Q::Acb,
    Q_dκ::Acb,
    Λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, σ) = Λ
    c = _c(κ, ϵ, Λ)

    # Compute enclosure of abs(Q)^2σ * Q differentiated w.r.t κ at ξ = ξ₁
    Q2σQ_dκ = let
        a = ArbSeries((real(Q), real(Q_dκ)))
        b = ArbSeries((imag(Q), imag(Q_dκ)))

        Q2σQ = abspow(a^2 + b^2, σ) * (a + im * b)

        Q2σQ[1]
    end

    I_P_dκ_2_1 = exp(-c * ξ₁^2) * F.P * ξ₁^(d - 2) * Q2σQ_dκ

    # Compute bound of hat_I_P_dκ_2. This is based on Lemma
    # REF(lemma:I_P-remainder-bounds-1)
    α = (2σ + 1) * v - 2 / σ + d # The lemma uses α to denote this value
    @assert α - 2 < 0 # Requirement for lemma

    hat_I_P_dκ_2_2_bound =
        (
            (2σ + 1) * (C.P_dξ + abs(d - 2) * C.P) / abs(α - 4) *
            norms.Q^2σ *
            norms.Q_dκ *
            ξ₁^-1 +
            (2σ + 1) * C.P / abs(α - 3) *
            (2σ * norms.Q_dξ * norms.Q_dκ + norms.Q * norms.Q_dκ_dξ) *
            norms.Q^(2σ - 1)
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^(α - 3)

    main = B_W(κ, ϵ, Λ) * (I_P_dκ_2_1 / 2c)
    remainder = add_error(zero(γ), abs(B_W(κ, ϵ, Λ) / 2c) * hat_I_P_dκ_2_2_bound)

    return main + remainder
end

function I_P_dϵ_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Q::Acb,
    Q_dξ::Acb,
    Q_dϵ::Acb,
    Λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    return I_P_dϵ_1_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dξ, Λ, F, C, norms) +
           I_P_dϵ_2_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dϵ, Λ, F, C, norms)
end

"""
    I_P_dϵ_1_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dξ, Λ, F, C, norms)

Compute an enclosure of ``I_P_dϵ_1`` using the expansion from Lemma
REF(lemma:I_P-derivatives-expansion-2).
"""
function I_P_dϵ_1_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Q::Acb,
    Q_dξ::Acb,
    Λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, σ) = Λ
    c = _c(κ, ϵ, Λ)

    # Compute enclosure of abs(Q)^2σ * Q and its first derivative at ξ = ξ₁
    Q2σQ, Q2σQ_dξ = let
        a = ArbSeries((real(Q), real(Q_dξ)))
        b = ArbSeries((imag(Q), imag(Q_dξ)))

        Q2σQ = abspow(a^2 + b^2, σ) * (a + im * b)

        Q2σQ[0], Q2σQ[1]
    end

    # For computing I_P_dϵ_1_2 we explicitly expand the derivative in
    # the formulas in Lemma REF(lemma:I_P-derivatives-expansion-2)

    I_P_dϵ_1_1 = exp(-c * ξ₁^2) * F.H * ξ₁^d * Q2σQ

    I_P_dϵ_1_2 =
        exp(-c * ξ₁^2) * (
            F.H_dξ * ξ₁^(d - 1) * Q2σQ +
            d * F.H * ξ₁^(d - 2) * Q2σQ +
            F.H * ξ₁^(d - 1) * Q2σQ_dξ
        )

    # Compute bound of hat_I_P_dϵ_1_3. This is based on Lemma
    # REF(lemma:I_P-remainder-bounds-2)
    α = (2σ + 1) * v - 2 / σ + d # The lemma uses α to denote this value
    @assert α - 2 < 0 # Requirement for lemma

    hat_I_P_dϵ_1_3_bound =
        (
            (C.H_dξ_dξ + abs(2d - 1) * C.H_dξ + abs(d * (d - 2)) * C.H) / abs(α - 4) *
            norms.Q^(2σ + 1) *
            ξ₁^-2 +
            (2σ + 1) * (2C.H_dξ + abs(2d - 1) * C.H) / abs(α - 3) *
            norms.Q^2σ *
            norms.Q_dξ *
            ξ₁^-1 +
            (2σ + 1) * C.H / abs(α - 2) *
            (2σ * norms.Q_dξ^2 + norms.Q * norms.Q_dξ_dξ) *
            norms.Q^(2σ - 1)
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^(α - 2)

    main = I_P_dϵ_1_1 / 2c + I_P_dϵ_1_2 / (2c)^2
    remainder = add_error(zero(γ), abs(1 / (2c)^2) * hat_I_P_dϵ_1_3_bound)

    return main + remainder
end

"""
    I_P_dϵ_2_enclose(γ, κ, ϵ, ξ₁, v, Q, Q_dϵ, Λ, F, C, norms)

Compute an enclosure of ``I_P_dϵ_2`` using the expansion from Lemma
REF(lemma:I_P-derivatives-expansion-1).
"""
function I_P_dϵ_2_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Q::Acb,
    Q_dϵ::Acb,
    Λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, σ) = Λ
    c = _c(κ, ϵ, Λ)

    # Compute enclosure of abs(Q)^2σ * Q differentiated w.r.t ϵ at ξ = ξ₁
    Q2σQ_dϵ = let
        a = ArbSeries((real(Q), real(Q_dϵ)))
        b = ArbSeries((imag(Q), imag(Q_dϵ)))

        Q2σQ = abspow(a^2 + b^2, σ) * (a + im * b)

        Q2σQ[1]
    end

    I_P_dϵ_2_1 = exp(-c * ξ₁^2) * F.P * ξ₁^(d - 2) * Q2σQ_dϵ

    # Compute bound of hat_I_P_dϵ_2. This is based on Lemma
    # REF(lemma:I_P-remainder-bounds-1)
    α = (2σ + 1) * v - 2 / σ + d # The lemma uses α to denote this value
    @assert α - 2 < 0 # Requirement for lemma

    hat_I_P_dϵ_2_2_bound =
        (
            (2σ + 1) * (C.P_dξ + abs(d - 2) * C.P) / abs(α - 4) *
            norms.Q^2σ *
            norms.Q_dϵ *
            ξ₁^-1 +
            (2σ + 1) * C.P / abs(α - 3) *
            (2σ * norms.Q_dξ * norms.Q_dϵ + norms.Q * norms.Q_dϵ_dξ) *
            norms.Q^(2σ - 1)
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^(α - 3)

    main = B_W(κ, ϵ, Λ) * (I_P_dϵ_2_1 / 2c)
    remainder = add_error(zero(γ), abs(B_W(κ, ϵ, Λ) / 2c) * hat_I_P_dϵ_2_2_bound)

    return main + remainder
end
