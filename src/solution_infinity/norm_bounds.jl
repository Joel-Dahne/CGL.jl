function C_u_dξ(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_E_dξ(κ, λ, ξ₁) * C_I_P(κ, ξ₁, v, λ) +
           (
        C_P_dξ(κ, λ, ξ₁) * C_I_E(κ, ξ₁, v, λ) +
        C_P(κ, λ, ξ₁) * C_J_E(κ, ξ₁, λ) +
        C_E(κ, λ, ξ₁) * C_J_P(κ, ξ₁, λ)
    ) * ξ₁^(-2)
end

function C_u_dξ_dξ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_E_dξ_dξ(κ, λ, ξ₁) * C_I_P_1(κ, ξ₁, v, λ) +
           C_E_dξ(κ, λ, ξ₁) * C_I_P_dξ(κ, ξ₁, v, λ) +
           C_E_dξ(κ, λ, ξ₁) * C_J_P(κ, ξ₁, λ) +
           C_E(κ, λ, ξ₁) * C_J_P_dξ(κ, ξ₁, λ) +
           (
               C_P_dξ_dξ(κ, λ, ξ₁) * C_I_E(κ, ξ₁, v, λ) +
               C_P_dξ(κ, λ, ξ₁) * C_I_E_dξ(κ, ξ₁, v, λ) +
               C_P_dξ(κ, λ, ξ₁) * C_J_E(κ, ξ₁, λ) +
               C_P(κ, λ, ξ₁) * C_J_E_dξ(κ, ξ₁, λ) +
               C_E(κ, λ, ξ₁) * C_J_P(κ, ξ₁, λ)
           ) * ξ₁^(-2)
end

function C_u_dξ_dξ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_E_dξ_dξ(κ, λ, ξ₁) * C_I_P_2(κ, ξ₁, v, λ) +
           (2λ.σ + 1) *
           (C_P(κ, λ, ξ₁) * C_J_E(κ, ξ₁, λ) + C_E(κ, λ, ξ₁) * C_J_P(κ, ξ₁, λ)) *
           ξ₁^(-2)
end

function C_u_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_P_dκ(κ, λ, ξ₁) * exp(-one(κ)) / v
end

function C_u_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    (; σ) = λ

    return C_P_dκ(κ, λ, ξ₁) * C_I_E(κ, ξ₁, v, λ) * exp(-one(κ)) / v * ξ₁^((2σ + 1)v - 2) +
           C_P(κ, λ, ξ₁) * C_I_E_dκ_1(κ, ξ₁, v, λ) * ξ₁^(2σ * v - 2) +
           C_P(κ, λ, ξ₁) * C_I_E_dκ_3(κ, ξ₁, v, λ) * log(ξ₁) * ξ₁^(2σ * v - 2) +
           C_E_dκ(κ, λ, ξ₁) * C_I_P_1(κ, ξ₁, v, λ) * ξ₁^(2σ * v - 2) +
           C_E(κ, λ, ξ₁) *
           (C_I_P_dκ_1(κ, ξ₁, v, λ) + C_I_P_dκ_2(κ, ξ₁, v, λ)) *
           ξ₁^(2σ * v - 2)
end

function C_u_dκ_3(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    (; σ) = λ

    return C_E_dκ(κ, λ, ξ₁) * C_I_P_2(κ, ξ₁, v, λ) * ξ₁^(2σ * v - 1) +
           C_E(κ, λ, ξ₁) *
           (C_I_P_dκ_3(κ, ξ₁, v, λ) + C_I_P_dκ_3(κ, ξ₁, v, λ) * ξ₁^(-1)) *
           ξ₁^(2σ * v - 2)
end

function C_u_dκ_4(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    (; σ) = λ

    return C_E(κ, λ, ξ₁) * C_I_P_dκ_5(κ, ξ₁, v, λ) * ξ₁^(2σ * v - 2)
end

function C_u_dκ_5(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    (; σ) = λ

    return C_E(κ, λ, ξ₁) * C_I_P_dκ_6(κ, ξ₁, v, λ) * ξ₁^(2σ * v - 2)
end

function C_u_dκ_6(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    (; σ) = λ

    return (
        C_P(κ, λ, ξ₁) * C_I_E_dκ_2(κ, ξ₁, v, λ) + C_E(κ, λ, ξ₁) * C_I_P_dκ_3(κ, ξ₁, v, λ)
    ) * ξ₁^(2σ * v - 2)
end
