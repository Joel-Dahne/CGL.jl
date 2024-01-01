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
    C::FunctionBounds,
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
            C.P_dξ_dξ_dξ / den1 * norm_u2σu * ξ₁^-3 +
            abs(3d - 9) * C.P_dξ_dξ / den1 * norm_u2σu * ξ₁^-3 +
            3C.P_dξ_dξ / den2 * norm_u2σu_dξ * ξ₁^-2 +
            abs(3d^2 - 21d + 33) * C.P_dξ / den1 * norm_u2σu * ξ₁^-3 +
            abs(6d - 18) * C.P_dξ / den2 * norm_u2σu_dξ * ξ₁^-2 +
            3C.P_dξ / den3 * norm_u2σu_dξ_dξ * ξ₁^-1 +
            abs((d - 2) * (d - 4) * (d - 6)) * C.P / den1 * norm_u2σu * ξ₁^-3 +
            abs(3d^2 - 21d + 33) * C.P / den2 * norm_u2σu_dξ * ξ₁^-2 +
            abs(3d - 9) * C.P / den3 * norm_u2σu_dξ_dξ * ξ₁^-1 +
            C.P / den4 * norm_u2σu_dξ_dξ_dξ
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
    norm_u_dξ::Arb,
    norm_u_dγ::Arb,
    norm_u_dξ_dγ::Arb,
    u::Acb,
    u_dγ::Acb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # Compute abs(u)^2σ * u differentiated w.r.t γ
    u2σu_dγ = let
        a = ArbSeries((real(u), real(u_dγ)))
        b = ArbSeries((imag(u), imag(u_dγ)))

        u2σu = abspow(a^2 + b^2, σ) * (a + im * b)

        u2σu[1]
    end

    I_P_dγ_1 = exp(-c * ξ₁^2) * P(ξ₁, (λ, κ)) * ξ₁^(d - 2) * u2σu_dγ

    # Compute bound of hat_I_P_dγ_2

    # Norms of required derivatives of abs(u)^2σ * u
    norm_u2σu_dγ = (2σ + 1) * norm_u^2σ * norm_u_dγ
    norm_u2σu_dξ_dγ =
        (2σ + 1) * norm_u^(2σ - 1) * (2σ * norm_u_dξ * norm_u_dγ + norm_u * norm_u_dξ_dγ)

    # Denominators coming from the integration
    den1 = abs((2σ + 1) * v - 2 / σ + d - 4)
    den2 = abs((2σ + 1) * v - 2 / σ + d - 3)

    hat_I_P_dγ_2_bound =
        (
            C.P_dξ / den1 * norm_u2σu_dγ * ξ₁^-1 +
            abs(d - 2) * C.P / den1 * norm_u2σu_dγ * ξ₁^-1 +
            C.P / den2 * norm_u2σu_dξ_dγ
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    main = B_W(κ, λ) * (I_P_dγ_1 / 2c)
    remainder = add_error(zero(γ), abs(B_W(κ, λ) / 2c) * hat_I_P_dγ_2_bound)

    return main + remainder
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
    u::Acb,
    u_dξ::Acb,
    u_dκ::Acb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    return I_P_dκ_1_enclose(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, u, u_dξ, λ, C) +
           I_P_dκ_2_enclose(
        γ,
        κ,
        ξ₁,
        v,
        norm_u,
        norm_u_dξ,
        norm_u_dκ,
        norm_u_dξ_dκ,
        u,
        u_dκ,
        λ,
        C,
    )
end

function I_P_dκ_1_enclose(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    u::Acb,
    u_dξ::Acb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # Compute abs(u)^2σ * u and its first derivative
    u2σu, u2σu_dξ = let
        a = ArbSeries((real(u), real(u_dξ)))
        b = ArbSeries((imag(u), imag(u_dξ)))

        u2σu = abspow(a^2 + b^2, σ) * (a + im * b)

        u2σu[0], u2σu[1]
    end

    I_P_dκ_1_1 = exp(-c * ξ₁^2) * D(ξ₁, (λ, κ)) * ξ₁^d * u2σu

    I_P_dκ_1_2 =
        exp(-c * ξ₁^2) * (
            D_dξ(ξ₁, (λ, κ)) * ξ₁^(d - 1) * u2σu +
            d * D(ξ₁, (λ, κ)) * ξ₁^(d - 2) * u2σu +
            D(ξ₁, (λ, κ)) * ξ₁^(d - 1) * u2σu_dξ
        )

    # Compute bound of hat_I_P_dκ_1_2

    # Norms of required derivatives of abs(u)^2σ * u
    norm_u2σu = norm_u^(2σ + 1)
    norm_u2σu_dξ = (2σ + 1) * norm_u^2σ * norm_u_dξ
    norm_u2σu_dξ_dξ =
        (2σ + 1) * norm_u^(2σ - 1) * (2σ * norm_u_dξ^2 + norm_u * norm_u_dξ_dξ)

    # Denominators coming from the integration
    den1 = abs((2σ + 1) * v - 2 / σ + d - 4)
    den2 = abs((2σ + 1) * v - 2 / σ + d - 3)
    den3 = abs((2σ + 1) * v - 2 / σ + d - 2)

    hat_I_P_dκ_1_3_bound =
        (
            C.D_dξ_dξ / den1 * norm_u2σu * ξ₁^-2 +
            abs(2d - 1) * C.D_dξ / den1 * norm_u2σu * ξ₁^-2 +
            2C.D_dξ / den2 * norm_u2σu_dξ * ξ₁^-1 +
            abs(d * (d - 2)) * C.D / den1 * norm_u2σu * ξ₁^-2 +
            abs(2d - 1) * C.D / den2 * norm_u2σu_dξ * ξ₁^-1 +
            C.D / den3 * norm_u2σu_dξ_dξ
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    main = I_P_dκ_1_1 / 2c + I_P_dκ_1_2 / (2c)^2
    remainder = add_error(zero(γ), abs(1 / (2c)^2) * hat_I_P_dκ_1_3_bound)

    return main + remainder
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
    u::Acb,
    u_dκ::Acb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # Compute abs(u)^2σ * u differentiated w.r.t κ
    u2σu_dκ = let
        a = ArbSeries((real(u), real(u_dκ)))
        b = ArbSeries((imag(u), imag(u_dκ)))

        u2σu = abspow(a^2 + b^2, σ) * (a + im * b)

        u2σu[1]
    end

    I_P_dκ_2_1 = exp(-c * ξ₁^2) * P(ξ₁, (λ, κ)) * ξ₁^(d - 2) * u2σu_dκ

    # Compute bound of hat_I_P_dκ_2_2

    # Norms of required derivatives of abs(u)^2σ * u
    norm_u2σu_dκ = (2σ + 1) * norm_u^2σ * norm_u_dκ
    norm_u2σu_dξ_dκ =
        (2σ + 1) * norm_u^(2σ - 1) * (2σ * norm_u_dξ * norm_u_dκ + norm_u * norm_u_dξ_dκ)

    # Denominators coming from the integration
    den1 = abs((2σ + 1) * v - 2 / σ + d - 4)
    den2 = abs((2σ + 1) * v - 2 / σ + d - 3)

    hat_I_P_dκ_2_2_bound =
        (
            C.P_dξ / den1 * norm_u2σu_dκ * ξ₁^-1 +
            abs(d - 2) * C.P / den1 * norm_u2σu_dκ * ξ₁^-1 +
            C.P / den2 * norm_u2σu_dξ_dκ
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    main = B_W(κ, λ) * (I_P_dκ_2_1 / 2c)
    remainder = add_error(zero(γ), abs(B_W(κ, λ) / 2c) * hat_I_P_dκ_2_2_bound)

    return main + remainder
end
