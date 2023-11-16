function I_P_1(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    u::Acb,
    u_dξ::Acb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # IMPROVE: The name of the function clashes with I_P_1 here. It
    # works, but looks awkward.

    I_P_1 = exp(-c * ξ₁^2) * P(ξ₁, (λ, κ)) * ξ₁^(d - 2) * abs(u)^2σ * u

    I_P_22 =
        exp(-c * ξ₁^2) * (
            P_dξ(ξ₁, (λ, κ)) * ξ₁^(d - 3) * abs(u)^2σ * u +
            (d - 2) * P(ξ₁, (λ, κ)) * ξ₁^(d - 4) * abs(u)^2σ * u +
            P(ξ₁, (λ, κ)) *
            ξ₁^(d - 3) *
            abs(u)^(2σ - 2) *
            (2σ * real(conj(u) * u_dξ) * u + abs(u)^2 * u_dξ)
        )

    # TODO: Improve the below bounds

    # These are like the non-alt versions but removing the factor
    # B_W(κ, λ) / 2c and the terms related to I_P_22
    C_I_P_2_2_alt =
        (
            C_P_dξ_dξ(κ, λ, ξ₁) +
            abs(2d - 5) * C_P_dξ(κ, λ, ξ₁) +
            abs((d - 2) * (d - 4)) * C_P(κ, λ, ξ₁)
        ) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 6)
    C_I_P_2_3_alt =
        (2σ + 1) * (
            (2C_P_dξ(κ, λ, ξ₁) + abs(2d - 5) * C_P(κ, λ, ξ₁)) /
            abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 5)
        )
    C_I_P_2_4_alt = (2σ + 1) * 2σ * C_P(κ, λ, ξ₁) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 4)
    C_I_P_2_5_alt = (2σ + 1) * C_P(κ, λ, ξ₁) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 4)
    bound_3 =
        (
            C_I_P_2_2_alt * norm_u^2 * ξ₁^-2 +
            C_I_P_2_3_alt * norm_u * norm_u_dξ * ξ₁^-1 +
            C_I_P_2_4_alt * norm_u_dξ^2 +
            C_I_P_2_5_alt * norm_u_dξ_dξ
        ) *
        norm_u^(2σ - 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 4)

    I_P_3 = add_error(zero(γ), bound_3)

    return B_W(κ, λ) / 2c * (I_P_1 + (I_P_22 + I_P_3) / 2c)
end
