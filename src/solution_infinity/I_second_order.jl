function R_u_sigma(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    R_u_bound_1 = zero(κ) # I_E is identically zero at ξ₁
    R_u_bound_2 =
        C_E(κ, λ, ξ₁) *
        (C_I_P_1_1(κ, ξ₁, v, λ) * norm_u * ξ₁^(-1) + C_I_P_1_2(κ, ξ₁, v, λ) * norm_u_dξ) *
        norm_u^2σ *
        ξ₁^((2σ + 1) * v + d - 4)

    R_u_bound = R_u_bound_1 + R_u_bound_2

    p0 = p_P(0, κ, λ)
    C_P_bound = C_P(1, κ, λ, ξ₁)

    abs(p0) > C_P_bound * ξ₁^(-2) || error("requires abs(p0) > R_P * ξ₁^(-2)")

    bound =
        R_u_bound * (
            norm_u^2σ +
            2σ *
            abs(γ)^2σ *
            C_P(κ, λ, ξ₁)^2σ *
            ξ₁^(-2σ * v) *
            (1 + 2σ * R_u_bound / (abs(γ) * (abs(p0) - C_P_bound * ξ₁^(-2))))
        )

    return add_error(zero(γ), bound)
end

function I_P_R(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; σ, d) = λ
    _, _, c = _abc(κ, λ)

    p0 = p_P(0, κ, λ)

    C_P_sigma =
        C_P(1, κ, λ, ξ₁) * (
            abs(p0)^2σ * σ +
            abs(p0)^(2σ - 1) * σ^2 * C_P(κ, λ, ξ₁) * ξ₁^-2 +
            C_P(κ, λ, ξ₁) * ξ₁^-1
        )

    R_I_P_1_bound =
        C_P_sigma / abs(-2 / σ + d - 4) * exp(-real(c) * ξ₁^2) * ξ₁^(-2 / σ + d - 4)
    R_I_P_2_bound =
        C_P_sigma^2 / abs(-2 / σ + d - 6) * exp(-real(c) * ξ₁^2) * ξ₁^(-2 / σ + d - 6)

    R_I_P_3_bound = let
        F1 = (norm_u^2σ + 2σ * abs(γ)^2σ * C_P(κ, λ, ξ₁)^2σ * ξ₁^(-2σ * v))
        F2 =
            (4σ^2 * abs(γ)^(2σ - 1) * C_P(κ, λ, ξ₁) * ξ₁^(-2σ * v)) /
            (abs(p0) - C_P(1, κ, λ, ξ₁) * ξ₁^-2)

        C1 = C_P(κ, λ, ξ₁) * C_I_E(κ, ξ₁, v, λ) * norm_u^(2σ + 1)
        C2 =
            C_E(κ, λ, ξ₁) *
            (C_I_P_1_1(κ, ξ₁, v, λ) * norm_u * ξ₁^-1 + C_I_P_1_2(κ, ξ₁, v, λ) * norm_u_dξ) *
            norm_u^2σ

        @assert 2σ * v - 2 / σ + d - 2 < 0
        @assert (4σ + 1) * v - 2 / σ + d - 4 < 0
        @assert (4σ + 1) * v - 2 / σ + d - 4 < 0
        R31 =
            C1 *
            (
                ξ₁^((2σ + 1) * v - 2) / abs(2σ * v - 2 / σ + d - 2) -
                ξ₁^((2σ + 1) * v - 2) / abs((4σ + 1) * v - 2 / σ + d - 4)
            ) *
            ξ₁^(2σ * v - 2 / σ + d - 2) +
            C2 * ξ₁^((4σ + 1) * v - 2 / σ + 2d - 6) / abs((4σ + 1) * v - 2 / σ + 2d - 6)
        @assert 2σ * v - 2 / σ + d - 2 < 0
        @assert (4σ + 1) * v - 2 / σ + d - 4 < 0
        @assert (6σ + 2) * v - 2 / σ + 2d - 8 < 0
        @assert (4σ + 1) * v - 2 / σ + d - 6 < 0
        @assert (6σ + 2) * v - 2 / σ + 2d - 8 < 0
        @assert (6σ + 2) * v - 2 / σ + 3d - 10 < 0
        R32 =
            C1^2 *
            (
                ξ₁^(2(2σ + 1) * v - 4) / abs(2σ * v - 2 / σ + d - 2) -
                2ξ₁^(2(2σ + 1) * v - 4) / abs((4σ + 1) * v - 2 / σ + d - 4) +
                ξ₁^(2(2σ + 1) * v - 4) / abs((6σ + 2) * v - 2 / σ + 2d - 8)
            ) *
            ξ₁^(2σ * v - 2 / σ + d - 2) +
            C1 *
            C2 *
            (
                ξ₁^((2σ + 1) * v - 2) / abs((4σ + 1) * v - 2 / σ + d - 6) -
                ξ₁^((2σ + 1) * v - 2) / abs((6σ + 2) * v - 2 / σ + 2d - 8)
            ) *
            ξ₁^((4σ + 1) * v - 2 / σ + d - 6) +
            C2^2 * ξ₁^((6σ + 2) * v - 2 / σ + 3d - 10) / abs((6σ + 2) * v - 2 / σ + 3d - 10)

        C_J_P(κ, ξ₁, λ) * exp(-real(c) * ξ₁^2) * (F1 * R31 + F2 * R32)
    end

    bound =
        abs(γ)^(2σ + 1) *
        abs(B_W(κ, λ)) *
        (2abs(p0)^(σ + 1) * R_I_P_1_bound + R_I_P_2_bound) + R_I_P_3_bound

    return add_error(zero(γ), bound)
end

function I_E_dξ_R(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    λ::AbstractGLParams{Arb},
)
    return C_J_E(κ, ξ₁, λ) *
           R_u_sigma(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ) *
           ξ₁^(2λ.σ * v - 3)
end

function I_P_dξ_R(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; σ, d) = λ
    _, _, c = _abc(κ, λ)

    return C_J_P(κ, ξ₁, λ) *
           R_u_sigma(γ, κ, ξ₁, v, norm_u, norm_u_dξ, λ) *
           exp(-real(c) * ξ₁^2) *
           ξ₁^(2σ * v - 2 / σ + d - 3)
end
