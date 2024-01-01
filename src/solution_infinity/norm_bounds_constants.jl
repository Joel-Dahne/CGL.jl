function C_u_dξ(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.E_dξ * C_I_P(κ, ξ₁, v, λ, C) +
           (
        C.P_dξ * C_I_E(κ, ξ₁, v, λ, C) +
        C.P * C_I_E_dξ(κ, ξ₁, v, λ, C) +
        C.E * C_I_P_dξ(κ, ξ₁, v, λ, C)
    ) * ξ₁^(-2)
end

function C_u_dξ_dξ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.E_dξ_dξ * C_I_P_1_1(κ, ξ₁, v, λ, C) +
           2C.E_dξ * C_I_P_dξ(κ, ξ₁, v, λ, C) +
           C.E * C.J_P_dξ +
           (
               C.P_dξ_dξ * C_I_E(κ, ξ₁, v, λ, C) +
               2C.P_dξ * C_I_E_dξ(κ, ξ₁, v, λ, C) +
               C.P * C.J_E_dξ
           ) * ξ₁^(-2)
end

function C_u_dξ_dξ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.E_dξ_dξ * C_I_P_1_2(κ, ξ₁, v, λ, C) +
           (2λ.σ + 1) * (C.P * C.J_E + C.E * C.J_P) * ξ₁^(-2)
end

function C_u_dξ_dξ_dξ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.E_dξ_dξ_dξ * (C_I_P_2_1(κ, ξ₁, v, λ, C) + C_I_P_2_2(κ, ξ₁, v, λ, C) * ξ₁^-2) +
           3C.E_dξ_dξ * C_I_P_dξ(κ, ξ₁, v, λ, C) +
           3C.E_dξ * C.J_P_dξ +
           C.E * C.J_P_dξ_dξ +
           (
               C.P_dξ_dξ_dξ * C_I_E(κ, ξ₁, v, λ, C) +
               3C.P_dξ_dξ * C_I_E_dξ(κ, ξ₁, v, λ, C) +
               3C.P_dξ * C.J_E_dξ +
               C.P * C.J_E_dξ_dξ
           ) * ξ₁^(-4)
end

function C_u_dξ_dξ_dξ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.E_dξ_dξ_dξ * C_I_P_2_3(κ, ξ₁, v, λ, C) +
           (2λ.σ + 1) * (
        3C.E_dξ * C.J_P + 2C.E * C.J_P_dξ + (3C.P_dξ * C.J_E + 2C.P * C.J_E_dξ) * ξ₁^(-2)
    )
end

function C_u_dξ_dξ_dξ_3(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.E_dξ_dξ_dξ * C_I_P_2_4(κ, ξ₁, v, λ, C) +
           (2λ.σ + 1) * 2λ.σ * (C.E * C.J_P + C.P * C.J_E) * ξ₁^-2
end

function C_u_dξ_dξ_dξ_4(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.E_dξ_dξ_dξ * C_I_P_2_4(κ, ξ₁, v, λ, C) +
           (2λ.σ + 1) * (C.E * C.J_P + C.P * C.J_E) * ξ₁^-2
end

function C_u_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.P_dκ * exp(-one(κ)) / v
end

function C_u_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = λ

    return C.P_dκ * C_I_E(κ, ξ₁, v, λ, C) * exp(-one(κ)) / v * ξ₁^((2σ + 1)v - 2) +
           C.P * C_I_E_dκ_1(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 2) +
           C.P * C_I_E_dκ_3(κ, ξ₁, v, λ, C) * log(ξ₁) * ξ₁^(2σ * v - 2) +
           C.E_dκ * C_I_P_1_1(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 2) +
           C.E *
           (
               C_I_P_dκ_1(κ, ξ₁, v, λ, C) +
               C_I_P_dκ_2(κ, ξ₁, v, λ, C) * ξ₁^-2 +
               C_I_P_dκ_3(κ, ξ₁, v, λ, C) * ξ₁^-2
           ) *
           ξ₁^(2σ * v - 2)
end

function C_u_dκ_3(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = λ

    return C.E_dκ * C_I_P_1_2(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 1) +
           C.E * (C_I_P_dκ_4(κ, ξ₁, v, λ, C) + C_I_P_dκ_5(κ, ξ₁, v, λ, C)) * ξ₁^(2σ * v - 3)
end

function C_u_dκ_4(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = λ

    return C.E * C_I_P_dκ_6(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 2)
end

function C_u_dκ_5(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = λ

    return C.E * C_I_P_dκ_7(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 2)
end

function C_u_dκ_6(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = λ

    return (C.P * C_I_E_dκ_2(κ, ξ₁, v, λ, C) + C.E * C_I_P_dκ_2_1(κ, ξ₁, v, λ, C)) *
           ξ₁^(2σ * v - 2)
end

function C_u_dξ_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = λ

    return C.P_dξ_dκ * C_I_E(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 3) +
           C.P_dξ *
           (C_I_E_dκ_1(κ, ξ₁, v, λ, C) + C_I_E_dκ_3(κ, ξ₁, v, λ, C) * log(ξ₁)) *
           ξ₁^(2σ * v - 3) +
           C.P_dξ * C.J_E * log(ξ₁) * ξ₁^(2σ * v - 3) +
           C.P * C.J_E_dκ * log(ξ₁) * ξ₁^(2σ * v - 3) +
           C.E_dξ_dκ *
           (C_I_P_2_1(κ, ξ₁, v, λ, C) + C_I_P_2_2(κ, ξ₁, v, λ, C) * ξ₁^-2) *
           ξ₁^(2σ * v - 1) +
           C.E_dξ *
           (C_I_P_dκ_1(κ, ξ₁, v, λ, C) + C_I_P_dκ_2(κ, ξ₁, v, λ, C)) *
           ξ₁^(2σ * v - 1) +
           C.E_dκ * C.J_P * ξ₁^(2σ * v - 1) +
           C.E * C.J_P_dκ * ξ₁^(2σ * v - 1)
end

function C_u_dξ_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = λ

    return C.P_dξ * C_I_E_dκ_2(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 3) +
           (2σ + 1) * C.P * C.J_E * ξ₁^(2σ * v - 3) +
           C.E_dξ * C_I_P_dκ_2_1(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 3) +
           (2σ + 1) * C.E * C.J_P * ξ₁^(2σ * v - 3)
end

function C_u_dξ_dκ_3(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = λ

    return C.E_dξ_dκ * C_I_P_2_3(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 2) +
           C.E_dξ *
           (C_I_P_dκ_3(κ, ξ₁, v, λ, C) + C_I_P_dκ_4(κ, ξ₁, v, λ, C) * ξ₁^-1) *
           ξ₁^(2σ * v - 1)
end

function C_u_dξ_dκ_4(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = λ

    return C.E_dξ_dκ * C_I_P_2_4(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 1) +
           C.E_dξ * C_I_P_dκ_5(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 1)
end

function C_u_dξ_dκ_5(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = λ

    return C.E_dξ_dκ * C_I_P_2_5(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 1) +
           C.E_dξ * C_I_P_dκ_6(κ, ξ₁, v, λ, C) * ξ₁^(2σ * v - 1)
end
