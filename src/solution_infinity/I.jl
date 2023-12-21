function I_P_enclose(
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
    λ::CGLParams{Arb},
)
    (; d, ω, σ, ϵ, δ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # Compute abs(u)^2σ * u and its first two derivatives
    u2σu, u2σu_dξ, u2σu_dξ_dξ = let
        # Enclose u_dξ_dξ using both the differential equation and the
        # bound for the norm and take the intersection.
        u_dξ_dξ_1 =
            -(
                (d - 1) / ξ₁ * u_dξ + im * κ * ξ₁ * u_dξ + im * κ / σ * u - ω * u +
                (1 + im * δ) * abs(u)^2σ * u
            ) / (1 - Acb(0, ϵ))

        u_dξ_dξ_2 = add_error(zero(γ), norm_u_dξ_dξ * ξ₁^(-1 / σ + v))

        u_dξ_dξ = Acb(
            Arblib.intersection(real(u_dξ_dξ_1), imag(u_dξ_dξ_2)),
            Arblib.intersection(imag(u_dξ_dξ_1), imag(u_dξ_dξ_2)),
        )

        a = ArbSeries((real(u), real(u_dξ), real(u_dξ_dξ)))
        b = ArbSeries((imag(u), imag(u_dξ), imag(u_dξ_dξ)))

        u2σu = abspow(a^2 + b^2, σ) * (a + im * b)

        u2σu[0], u2σu[1], u2σu[2]
    end

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

    main = B_W(κ, λ) * (I_P_1 / 2c + I_P_2 / (2c)^2 + I_P_3 / (2c)^3)
    remainder = add_error(zero(γ), abs(B_W(κ, λ) / (2c)^3) * hat_I_P_4_bound)

    return main + remainder
end

function I_P_dγ_enclose(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    λ::CGLParams{Arb},
)
    bound = I_P_dγ_bound_1(κ, ξ₁, v, norm_u, norm_u_dγ, λ)

    return add_error(zero(γ), bound)
end

function I_P_dγ_enclose(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dγ::Arb,
    λ::CGLParams{Arb},
)
    # The separate bounds are given in the paper. We compute both and
    # take the minimum.
    bound = min(
        I_P_dγ_bound_1(κ, ξ₁, v, norm_u, norm_u_dγ, λ),
        I_P_dγ_bound_2(κ, ξ₁, v, norm_u, norm_u_dγ, norm_u_dξ, norm_u_dξ_dγ, λ),
    )

    return add_error(zero(γ), bound)
end

function I_P_dγ_bound_1(
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    λ::CGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0 # Required for integral to converge

    bound =
        (2σ + 1) *
        C_I_P(κ, ξ₁, v, λ) *
        norm_u^2σ *
        norm_u_dγ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    return bound
end

# Similar to the above method but doing one step of partial integration
function I_P_dγ_bound_2(
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dγ::Arb,
    λ::CGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < 0 # Required for integral to converge

    bound =
        (
            (2σ + 1) * C_I_P_1_1(κ, ξ₁, v, λ) * norm_u * norm_u_dγ * ξ₁^(-1) +
            C_I_P_1_2(κ, ξ₁, v, λ) * (2σ * norm_u_dξ * norm_u_dγ + norm_u_dξ_dγ)
        ) *
        norm_u^(2σ - 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    return bound
end

function I_P_dκ_enclose(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    norm_u_dκ::Arb,
    norm_u_dξ_dκ::Arb,
    λ::CGLParams{Arb},
)
    return I_P_dκ_1_enclose(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dκ, λ) +
           I_P_dκ_2_enclose(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dκ, norm_u_dξ_dκ, λ)
end

function I_P_dκ_1_enclose(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    norm_u_dκ::Arb,
    λ::CGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    bound1 =
        C_I_P_dκ_1(κ, ξ₁, v, λ) *
        norm_u^(2σ + 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    bound2 =
        (C_I_P_dκ_2(κ, ξ₁, v, λ) * norm_u * ξ₁^-1 + C_I_P_dκ_4(κ, ξ₁, v, λ) * norm_u_dξ) *
        norm_u^2σ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    bound3 =
        (
            C_I_P_dκ_3(κ, ξ₁, v, λ) * norm_u^2 * ξ₁^-2 +
            C_I_P_dκ_5(κ, ξ₁, v, λ) * norm_u * norm_u_dξ * ξ₁^-1 +
            C_I_P_dκ_6(κ, ξ₁, v, λ) * norm_u_dξ^2 +
            C_I_P_dκ_7(κ, ξ₁, v, λ) * norm_u * norm_u_dξ_dξ
        ) *
        norm_u^(2σ - 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    bound = bound1 + bound2 + bound3

    return add_error(zero(γ), bound)
end

function I_P_dκ_2_enclose(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dκ::Arb,
    norm_u_dξ_dκ::Arb,
    λ::CGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    bound =
        (
            (2σ + 1) * C_I_P_1_1(κ, ξ₁, v, λ) * norm_u * norm_u_dκ * ξ₁^(-1) +
            C_I_P_1_2(κ, ξ₁, v, λ) * (2σ * norm_u_dξ * norm_u_dκ + norm_u_dξ_dκ)
        ) *
        norm_u^(2σ - 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    return add_error(zero(γ), bound)
end
