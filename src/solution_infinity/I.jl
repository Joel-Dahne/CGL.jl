function I_P_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    u::Acb,
    u_dξ::Acb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, ω, σ, δ) = λ

    _, _, c = _abc(κ, ϵ, λ)

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

        u_dξ_dξ_2 = add_error(zero(γ), norms.Q_dξ_dξ * ξ₁^(-1 / σ + v))

        u_dξ_dξ = Acb(
            Arblib.intersection(real(u_dξ_dξ_1), imag(u_dξ_dξ_2)),
            Arblib.intersection(imag(u_dξ_dξ_1), imag(u_dξ_dξ_2)),
        )

        a = ArbSeries((real(u), real(u_dξ), real(u_dξ_dξ)))
        b = ArbSeries((imag(u), imag(u_dξ), imag(u_dξ_dξ)))

        u2σu = abspow(a^2 + b^2, σ) * (a + im * b)

        u2σu[0], u2σu[1], 2u2σu[2]
    end

    I_P_1 = exp(-c * ξ₁^2) * F.P * ξ₁^(d - 2) * u2σu

    I_P_2 =
        exp(-c * ξ₁^2) * (
            F.P_dξ * ξ₁^(d - 3) * u2σu +
            (d - 2) * F.P * ξ₁^(d - 4) * u2σu +
            F.P * ξ₁^(d - 3) * u2σu_dξ
        )

    I_P_3 =
        exp(-c * ξ₁^2) * (
            F.P_dξ_dξ * ξ₁^(d - 4) * u2σu +
            (2d - 5) * F.P_dξ * ξ₁^(d - 5) * u2σu +
            2F.P_dξ * ξ₁^(d - 4) * u2σu_dξ +
            (d - 2) * (d - 4) * F.P * ξ₁^(d - 6) * u2σu +
            (2d - 5) * F.P * ξ₁^(d - 5) * u2σu_dξ +
            F.P * ξ₁^(d - 5) * u2σu_dξ_dξ
        )

    # Compute bound of hat_I_P_4

    # Denominators coming from the integration
    den1 = abs((2σ + 1) * v - 2 / σ + d - 8)
    den2 = abs((2σ + 1) * v - 2 / σ + d - 7)
    den3 = abs((2σ + 1) * v - 2 / σ + d - 6)
    den4 = abs((2σ + 1) * v - 2 / σ + d - 5)

    hat_I_P_4_bound =
        (
            C.P_dξ_dξ_dξ / den1 * norms.Q2σQ * ξ₁^-3 +
            abs(3d - 9) * C.P_dξ_dξ / den1 * norms.Q2σQ * ξ₁^-3 +
            3C.P_dξ_dξ / den2 * norms.Q2σQ_dξ * ξ₁^-2 +
            abs(3d^2 - 21d + 33) * C.P_dξ / den1 * norms.Q2σQ * ξ₁^-3 +
            abs(6d - 18) * C.P_dξ / den2 * norms.Q2σQ_dξ * ξ₁^-2 +
            3C.P_dξ / den3 * norms.Q2σQ_dξ_dξ * ξ₁^-1 +
            abs((d - 2) * (d - 4) * (d - 6)) * C.P / den1 * norms.Q2σQ * ξ₁^-3 +
            abs(3d^2 - 21d + 33) * C.P / den2 * norms.Q2σQ_dξ * ξ₁^-2 +
            abs(3d - 9) * C.P / den3 * norms.Q2σQ_dξ_dξ * ξ₁^-1 +
            C.P / den4 * norms.Q2σQ_dξ_dξ_dξ
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 5)

    main = B_W(κ, ϵ, λ) * (I_P_1 / 2c + I_P_2 / (2c)^2 + I_P_3 / (2c)^3)
    remainder = add_error(zero(γ), abs(B_W(κ, ϵ, λ) / (2c)^3) * hat_I_P_4_bound)

    return main + remainder
end

function I_P_dγ_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    u::Acb,
    u_dγ::Acb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, σ) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # Compute abs(u)^2σ * u differentiated w.r.t γ
    u2σu_dγ = let
        a = ArbSeries((real(u), real(u_dγ)))
        b = ArbSeries((imag(u), imag(u_dγ)))

        u2σu = abspow(a^2 + b^2, σ) * (a + im * b)

        u2σu[1]
    end

    I_P_dγ_1 = exp(-c * ξ₁^2) * F.P * ξ₁^(d - 2) * u2σu_dγ

    # Compute bound of hat_I_P_dγ_2

    # Denominators coming from the integration
    den1 = abs((2σ + 1) * v - 2 / σ + d - 4)
    den2 = abs((2σ + 1) * v - 2 / σ + d - 3)

    hat_I_P_dγ_2_bound =
        (
            C.P_dξ / den1 * norms.Q2σQ_dγ * ξ₁^-1 +
            abs(d - 2) * C.P / den1 * norms.Q2σQ_dγ * ξ₁^-1 +
            C.P / den2 * norms.Q2σQ_dγ_dξ
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    main = B_W(κ, ϵ, λ) * (I_P_dγ_1 / 2c)
    remainder = add_error(zero(γ), abs(B_W(κ, ϵ, λ) / 2c) * hat_I_P_dγ_2_bound)

    return main + remainder
end

function I_P_dκ_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    u::Acb,
    u_dξ::Acb,
    u_dκ::Acb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    return I_P_dκ_1_enclose(γ, κ, ϵ, ξ₁, v, u, u_dξ, λ, F, C, norms) +
           I_P_dκ_2_enclose(γ, κ, ϵ, ξ₁, v, u, u_dκ, λ, F, C, norms)
end

function I_P_dκ_1_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    u::Acb,
    u_dξ::Acb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, σ) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # Compute abs(u)^2σ * u and its first derivative
    u2σu, u2σu_dξ = let
        a = ArbSeries((real(u), real(u_dξ)))
        b = ArbSeries((imag(u), imag(u_dξ)))

        u2σu = abspow(a^2 + b^2, σ) * (a + im * b)

        u2σu[0], u2σu[1]
    end

    I_P_dκ_1_1 = exp(-c * ξ₁^2) * F.D * ξ₁^d * u2σu

    I_P_dκ_1_2 =
        exp(-c * ξ₁^2) * (
            F.D_dξ * ξ₁^(d - 1) * u2σu +
            d * F.D * ξ₁^(d - 2) * u2σu +
            F.D * ξ₁^(d - 1) * u2σu_dξ
        )

    # Compute bound of hat_I_P_dκ_1_2

    # Denominators coming from the integration
    den1 = abs((2σ + 1) * v - 2 / σ + d - 4)
    den2 = abs((2σ + 1) * v - 2 / σ + d - 3)
    den3 = abs((2σ + 1) * v - 2 / σ + d - 2)

    hat_I_P_dκ_1_3_bound =
        (
            C.D_dξ_dξ / den1 * norms.Q2σQ * ξ₁^-2 +
            abs(2d - 1) * C.D_dξ / den1 * norms.Q2σQ * ξ₁^-2 +
            2C.D_dξ / den2 * norms.Q2σQ_dξ * ξ₁^-1 +
            abs(d * (d - 2)) * C.D / den1 * norms.Q2σQ * ξ₁^-2 +
            abs(2d - 1) * C.D / den2 * norms.Q2σQ_dξ * ξ₁^-1 +
            C.D / den3 * norms.Q2σQ_dξ_dξ
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
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    u::Acb,
    u_dκ::Acb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, σ) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # Compute abs(u)^2σ * u differentiated w.r.t κ
    u2σu_dκ = let
        a = ArbSeries((real(u), real(u_dκ)))
        b = ArbSeries((imag(u), imag(u_dκ)))

        u2σu = abspow(a^2 + b^2, σ) * (a + im * b)

        u2σu[1]
    end

    I_P_dκ_2_1 = exp(-c * ξ₁^2) * F.P * ξ₁^(d - 2) * u2σu_dκ

    # Compute bound of hat_I_P_dκ_2_2

    # Denominators coming from the integration
    den1 = abs((2σ + 1) * v - 2 / σ + d - 4)
    den2 = abs((2σ + 1) * v - 2 / σ + d - 3)

    hat_I_P_dκ_2_2_bound =
        (
            C.P_dξ / den1 * norms.Q2σQ_dκ * ξ₁^-1 +
            abs(d - 2) * C.P / den1 * norms.Q2σQ_dκ * ξ₁^-1 +
            C.P / den2 * norms.Q2σQ_dκ_dξ
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    main = B_W(κ, ϵ, λ) * (I_P_dκ_2_1 / 2c)
    remainder = add_error(zero(γ), abs(B_W(κ, ϵ, λ) / 2c) * hat_I_P_dκ_2_2_bound)

    return main + remainder
end

function I_P_dϵ_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    u::Acb,
    u_dξ::Acb,
    u_dϵ::Acb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    return I_P_dϵ_1_enclose(γ, κ, ϵ, ξ₁, v, u, u_dξ, λ, F, C, norms) +
           I_P_dϵ_2_enclose(γ, κ, ϵ, ξ₁, v, u, u_dϵ, λ, F, C, norms)
end

function I_P_dϵ_1_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    u::Acb,
    u_dξ::Acb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, σ) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # Compute abs(u)^2σ * u and its first derivative
    u2σu, u2σu_dξ = let
        a = ArbSeries((real(u), real(u_dξ)))
        b = ArbSeries((imag(u), imag(u_dξ)))

        u2σu = abspow(a^2 + b^2, σ) * (a + im * b)

        u2σu[0], u2σu[1]
    end

    I_P_dϵ_1_1 = exp(-c * ξ₁^2) * F.H * ξ₁^d * u2σu

    I_P_dϵ_1_2 =
        exp(-c * ξ₁^2) * (
            F.H_dξ * ξ₁^(d - 1) * u2σu +
            d * F.H * ξ₁^(d - 2) * u2σu +
            F.H * ξ₁^(d - 1) * u2σu_dξ
        )

    # Compute bound of hat_I_P_dϵ_1_2

    # Denominators coming from the integration
    den1 = abs((2σ + 1) * v - 2 / σ + d - 4)
    den2 = abs((2σ + 1) * v - 2 / σ + d - 3)
    den3 = abs((2σ + 1) * v - 2 / σ + d - 2)

    hat_I_P_dϵ_1_3_bound =
        (
            C.H_dξ_dξ / den1 * norms.Q2σQ * ξ₁^-2 +
            abs(2d - 1) * C.H_dξ / den1 * norms.Q2σQ * ξ₁^-2 +
            2C.H_dξ / den2 * norms.Q2σQ_dξ * ξ₁^-1 +
            abs(d * (d - 2)) * C.H / den1 * norms.Q2σQ * ξ₁^-2 +
            abs(2d - 1) * C.H / den2 * norms.Q2σQ_dξ * ξ₁^-1 +
            C.H / den3 * norms.Q2σQ_dξ_dξ
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    main = I_P_dϵ_1_1 / 2c + I_P_dϵ_1_2 / (2c)^2
    remainder = add_error(zero(γ), abs(1 / (2c)^2) * hat_I_P_dϵ_1_3_bound)

    return main + remainder
end

function I_P_dϵ_2_enclose(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    u::Acb,
    u_dϵ::Acb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures,
    C::FunctionBounds,
    norms::NormBounds,
)
    (; d, σ) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    # Compute abs(u)^2σ * u differentiated w.r.t κ
    u2σu_dϵ = let
        a = ArbSeries((real(u), real(u_dϵ)))
        b = ArbSeries((imag(u), imag(u_dϵ)))

        u2σu = abspow(a^2 + b^2, σ) * (a + im * b)

        u2σu[1]
    end

    I_P_dϵ_2_1 = exp(-c * ξ₁^2) * F.P * ξ₁^(d - 2) * u2σu_dϵ

    # Compute bound of hat_I_P_dϵ_2_2

    # Denominators coming from the integration
    den1 = abs((2σ + 1) * v - 2 / σ + d - 4)
    den2 = abs((2σ + 1) * v - 2 / σ + d - 3)

    hat_I_P_dϵ_2_2_bound =
        (
            C.P_dξ / den1 * norms.Q2σQ_dϵ * ξ₁^-1 +
            abs(d - 2) * C.P / den1 * norms.Q2σQ_dϵ * ξ₁^-1 +
            C.P / den2 * norms.Q2σQ_dϵ_dξ
        ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    main = B_W(κ, ϵ, λ) * (I_P_dϵ_2_1 / 2c)
    remainder = add_error(zero(γ), abs(B_W(κ, ϵ, λ) / 2c) * hat_I_P_dϵ_2_2_bound)

    return main + remainder
end
