function I_P_1(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    norm_u_dξ_dξ_dξ::Arb,
    u::Acb,
    u_dξ::Acb,
    λ::AbstractGLParams{Arb},
)
    (; d, ω, σ, ϵ, δ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # Compute abs(u)^2σ * u and its first two derivatives
    u2σu, u2σu_dξ, u2σu_dξ_dξ = let
        # Use differential equation
        u_dξ_dξ =
            -(
                (d - 1) / ξ₁ * u_dξ + im * κ * ξ₁ * u_dξ + im * κ / σ * u - ω * u +
                (1 + im * δ) * abs(u)^2σ * u
            ) / (1 - Acb(0, ϵ))

        a = AcbSeries((real(u), real(u_dξ), real(u_dξ_dξ)))
        b = AcbSeries((imag(u), imag(u_dξ), imag(u_dξ_dξ)))

        u2σu = (a^2 + b^2)^σ * (a + im * b)

        u2σu[0], u2σu[1], u2σu[2]
    end

    # IMPROVE: The name of the function clashes with I_P_1 here. It
    # works, but looks awkward.

    I_P_1 = exp(-c * ξ₁^2) * P(ξ₁, (λ, κ)) * ξ₁^(d - 2) * u2σu

    I_P_2 =
        exp(-c * ξ₁^2) * (
            P_dξ(ξ₁, (λ, κ)) * ξ₁^(d - 3) * u2σu +
            (d - 2) * P(ξ₁, (λ, κ)) * ξ₁^(d - 4) * u2σu +
            P(ξ₁, (λ, κ)) * ξ₁^(d - 3) * u2σu_dξ
        )

    I_P_3 =
        exp(-c * ξ₁^2) * (
            P_dξ_dξ(ξ₁, (λ, κ)) * ξ₁^(d - 4) * u2σu +
            (2d - 5) * P_dξ(ξ₁, (λ, κ)) * ξ₁^(d - 5) * u2σu +
            2P_dξ(ξ₁, (λ, κ)) * ξ₁^(d - 4) * u2σu_dξ +
            (d - 2) * (d - 4) * P(ξ₁, (λ, κ)) * ξ₁^(d - 6) * u2σu +
            (2d - 5) * P(ξ₁, (λ, κ)) * ξ₁^(d - 5) * u2σu_dξ +
            P(ξ₁, (λ, κ)) * ξ₁^(d - 5) * u2σu_dξ_dξ
        )

    # Compute bound of hat_I_P_4

    # Norms of abs(u)^2σ * u and its derivatives
    norm_u2σu = norm_u^(2σ + 1)
    norm_u2σu_dξ = (2σ + 1) * norm_u^2σ * norm_u_dξ
    norm_u2σu_dξ_dξ =
        (2σ + 1) * norm_u^(2σ - 1) * (2σ * norm_u_dξ^2 + norm_u * norm_u_dξ_dξ)
    norm_u2σu_dξ_dξ_dξ =
        (2σ + 1) *
        norm_u^(2σ - 2) *
        (
            2σ * (2σ - 1) * norm_u_dξ^3 +
            (2σ + 2) * norm_u * norm_u_dξ * norm_u_dξ_dξ +
            norm_u^2 * norm_u_dξ_dξ_dξ
        )

    # Denominators coming from the integration
    den1 = abs((2σ + 1) * v - 2 / σ + d - 8)
    den2 = abs((2σ + 1) * v - 2 / σ + d - 7)
    den3 = abs((2σ + 1) * v - 2 / σ + d - 6)
    den4 = abs((2σ + 1) * v - 2 / σ + d - 5)

    hat_I_P_4_bound =
        (
            C_P_dξ_dξ_dξ(κ, λ, ξ₁) / den1 * norm_u2σu * ξ₁^-3 +
            (3d - 9) * C_P_dξ_dξ(κ, λ, ξ₁) / den1 * norm_u2σu * ξ₁^-3 +
            3C_P_dξ_dξ(κ, λ, ξ₁) / den2 * norm_u2σu_dξ * ξ₁^-2 +
            (3d^2 - 21d + 33) * C_P_dξ(κ, λ, ξ₁) / den1 * norm_u2σu * ξ₁^-3 +
            (6d - 18) * C_P_dξ(κ, λ, ξ₁) / den2 * norm_u2σu_dξ * ξ₁^-2 +
            3C_P_dξ(κ, λ, ξ₁) / den3 * norm_u2σu_dξ_dξ * ξ₁^-1 +
            (d - 2) * (d - 4) * (d - 6) * C_P(κ, λ, ξ₁) / den1 * norm_u2σu * ξ₁^-3 +
            (3d^2 - 21d + 33) * C_P(κ, λ, ξ₁) / den2 * norm_u2σu_dξ * ξ₁^-2 +
            (3d - 9) * C_P(κ, λ, ξ₁) / den3 * norm_u2σu_dξ_dξ * ξ₁^-1 +
            C_P(κ, λ, ξ₁) / den4 * norm_u2σu_dξ_dξ_dξ
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 5)

    hat_I_P_4 = add_error(zero(γ), hat_I_P_4_bound)

    return B_W(κ, λ) * (I_P_1 / 2c + I_P_2 / (2c)^2 + I_P_3 / (2c)^3 + hat_I_P_4 / (2c)^3)
end
