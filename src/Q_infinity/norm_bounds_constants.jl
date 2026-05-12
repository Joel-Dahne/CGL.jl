function C_T(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    return C.P * CI.I_E + C.E * CI.I_P
end

function C_Q_dξ(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    return C.E_dξ * CI.I_P + (C.P_dξ * CI.I_E + C.P * CI.I_E_dξ + C.E * CI.I_P_dξ) * ξ₁^-2
end

function C_Q_dξ_dξ_1(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    return C.E_dξ_dξ * CI.I_P_1_1 +
           2C.E_dξ * CI.I_P_dξ +
           C.E * C.J_P_dξ +
           (C.P_dξ_dξ * CI.I_E + 2C.P_dξ * CI.I_E_dξ + C.P * C.J_E_dξ) * ξ₁^-2
end

function C_Q_dξ_dξ_2(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return C.E_dξ_dξ * CI.I_P_1_2 + (2σ + 1) * (C.P * C.J_E + C.E * C.J_P) * ξ₁^-2
end

function C_Q_dξ_dξ_dξ_1(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    return C.E_dξ_dξ_dξ * CI.I_P_2_1 +
           3C.E_dξ_dξ * CI.I_P_dξ +
           3C.E_dξ * C.J_P_dξ +
           C.E * C.J_P_dξ_dξ +
           (
               C.P_dξ_dξ_dξ * CI.I_E +
               3C.P_dξ_dξ * CI.I_E_dξ +
               3C.P_dξ * C.J_E_dξ +
               C.P * C.J_E_dξ_dξ
           ) * ξ₁^-4
end

function C_Q_dξ_dξ_dξ_2(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return C.E_dξ_dξ_dξ * CI.I_P_2_2 +
           (2σ + 1) *
           (3C.E_dξ * C.J_P + 2C.E * C.J_P_dξ + (3C.P_dξ * C.J_E + 2C.P * C.J_E_dξ) * ξ₁^-2)
end

function C_Q_dξ_dξ_dξ_3(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return C.E_dξ_dξ_dξ * CI.I_P_2_3 + (2σ + 1) * 2σ * (C.E * C.J_P + C.P * C.J_E) * ξ₁^-2
end

function C_Q_dξ_dξ_dξ_4(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return C.E_dξ_dξ_dξ * CI.I_P_2_4 + (2σ + 1) * (C.E * C.J_P + C.P * C.J_E) * ξ₁^-2
end

function C_Q_dκ_1(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    return C.P_dκ * exp(-Arb(1)) / v
end

function C_Q_dκ_2(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (
        C.P_dκ * CI.I_E * exp(-Arb(1)) / v * ξ₁^v +
        C.P * CI.I_E_dκ +
        C.E_dκ * CI.I_P_1_1 +
        C.E * CI.I_P_dκ_1_1
    ) * ξ₁^(2σ * v - 2)
end

function C_Q_dκ_3(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (C.E_dκ * CI.I_P_1_2 + C.E * CI.I_P_dκ_1_2 * ξ₁^-2) * ξ₁^(2σ * v - 1)
end

function C_Q_dκ_4(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return C.E * CI.I_P_dκ_1_3 * ξ₁^(2σ * v - 2)
end

function C_Q_dκ_5(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return C.E * CI.I_P_dκ_1_4 * ξ₁^(2σ * v - 2)
end

function C_Q_dκ_6(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (2σ + 1) * (C.P * CI.I_E + C.E * CI.I_P) * ξ₁^(2σ * v - 2)
end

function C_Q_dξ_dκ_1(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (
        C.P_dξ_dκ * CI.I_E * log(ξ₁) * ξ₁^-2 +
        C.P_dκ * C.J_E * log(ξ₁) * ξ₁^-2 +
        C.P_dξ * CI.I_E_dκ * ξ₁^-2 +
        C.P * C.J_E_dκ * log(ξ₁) * ξ₁^-2 +
        C.E_dξ_dκ * CI.I_P_2_1 +
        C.E_dκ * C.J_P +
        C.E_dξ * CI.I_P_dκ_1_1 +
        C.E * C.J_P_dκ
    ) * ξ₁^(2σ * v - 1)
end

function C_Q_dξ_dκ_2(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (2σ + 1) *
           (C.P_dξ * CI.I_E + C.P * C.J_E + C.E_dξ * CI.I_P + C.E * C.J_P) *
           ξ₁^(2σ * v - 3)
end

function C_Q_dξ_dκ_3(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (C.E_dξ_dκ * CI.I_P_2_2 + C.E_dξ * CI.I_P_dκ_1_2) * ξ₁^(2σ * v - 2)
end

function C_Q_dξ_dκ_4(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (C.E_dξ_dκ * CI.I_P_2_3 + C.E_dξ * CI.I_P_dκ_1_3) * ξ₁^(2σ * v - 1)
end

function C_Q_dξ_dκ_5(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (C.E_dξ_dκ * CI.I_P_2_4 + C.E_dξ * CI.I_P_dκ_1_4) * ξ₁^(2σ * v - 1)
end

function C_Q_dϵ_1(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    return C.P_dϵ * ξ₁^-v
end

function C_Q_dϵ_2(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (C.P_dϵ * CI.I_E + C.P * CI.I_E_dϵ + C.E_dϵ * CI.I_P_1_1 + C.E * CI.I_P_dϵ_1_1) *
           ξ₁^(2σ * v - 2)
end

function C_Q_dϵ_3(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (C.E_dϵ * CI.I_P_1_2 + C.E * CI.I_P_dϵ_1_2 * ξ₁^-2) * ξ₁^(2σ * v - 1)
end

function C_Q_dϵ_4(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return C.E * CI.I_P_dϵ_1_3 * ξ₁^(2σ * v - 2)
end

function C_Q_dϵ_5(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return C.E * CI.I_P_dϵ_1_4 * ξ₁^(2σ * v - 2)
end

function C_Q_dϵ_6(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (2σ + 1) * (C.P * CI.I_E + C.E * CI.I_P) * ξ₁^(2σ * v - 2)
end

function C_Q_dξ_dϵ_1(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (
        C.P_dξ_dϵ * CI.I_E * ξ₁^-2 +
        C.P_dϵ * C.J_E * ξ₁^-2 +
        C.P_dξ * CI.I_E_dϵ * ξ₁^-2 +
        C.P * C.J_E_dϵ * ξ₁^-2 +
        C.E_dξ_dϵ * CI.I_P_2_1 +
        C.E_dϵ * C.J_P +
        C.E_dξ * CI.I_P_dϵ_1_1 +
        C.E * C.J_P_dϵ
    ) * ξ₁^(2σ * v - 1)
end

function C_Q_dξ_dϵ_2(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (2σ + 1) *
           (C.P_dξ * CI.I_E + C.P * C.J_E + C.E_dξ * CI.I_P + C.E * C.J_P) *
           ξ₁^(2σ * v - 3)
end

function C_Q_dξ_dϵ_3(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (C.E_dξ_dϵ * CI.I_P_2_2 + C.E_dξ * CI.I_P_dϵ_1_2) * ξ₁^(2σ * v - 2)
end

function C_Q_dξ_dϵ_4(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (C.E_dξ_dϵ * CI.I_P_2_3 + C.E_dξ * CI.I_P_dϵ_1_3) * ξ₁^(2σ * v - 1)
end

function C_Q_dξ_dϵ_5(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; σ) = Λ
    return (C.E_dξ_dϵ * CI.I_P_2_4 + C.E_dξ * CI.I_P_dϵ_1_4) * ξ₁^(2σ * v - 1)
end
