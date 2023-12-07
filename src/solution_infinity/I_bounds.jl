function C_I_E(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    @assert (2λ.σ + 1) * v - 2 < 0
    return C_J_E(κ, ξ₁, λ) / abs((2λ.σ + 1) * v - 2)
end

function C_I_P(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    @assert (2λ.σ + 1) * v - 2 / λ.σ + λ.d - 2 < 0
    return C_J_P(κ, ξ₁, λ) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 2)
end

function C_I_P_1_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 3 < 0

    _, _, c = _abc(κ, λ)

    return abs(B_W(κ, λ) / 2c) * (
        C_P(κ, λ, ξ₁) +
        C_P_dξ(κ, λ, ξ₁) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 3) +
        abs(d - 2) * C_P(κ, λ, ξ₁) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 4)
    )
end

function C_I_P_1_2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    @assert (2λ.σ + 1) * v - 2 / λ.σ + λ.d - 3 < 0

    (; σ, d) = λ

    _, _, c = _abc(κ, λ)

    return abs(B_W(κ, λ) / 2c) * (2σ + 1) * C_P(κ, λ, ξ₁) /
           abs((2σ + 1) * v - 2 / σ + d - 3)
end

function C_I_P_2_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    _, _, c = _abc(κ, λ)

    return abs(B_W(κ, λ) / 2c) * C_P(κ, λ, ξ₁)
end

function C_I_P_2_2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    _, _, c = _abc(κ, λ)

    return abs(B_W(κ, λ) / 4c^2) * (
        C_P_dξ(κ, λ, ξ₁) +
        abs(d - 2) * C_P(κ, λ, ξ₁) +
        (
            C_P_dξ_dξ(κ, λ, ξ₁) +
            abs(2d - 5) * C_P_dξ(κ, λ, ξ₁) +
            abs((d - 2) * (d - 4)) * C_P(κ, λ, ξ₁)
        ) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 6)
    )
end

function C_I_P_2_3(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    _, _, c = _abc(κ, λ)

    return abs(B_W(κ, λ) / 4c^2) *
           (2σ + 1) *
           (
               C_P(κ, λ, ξ₁) +
               (2C_P_dξ(κ, λ, ξ₁) + abs(2d - 5) * C_P(κ, λ, ξ₁)) /
               abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 5)
           )
end

function C_I_P_2_4(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    _, _, c = _abc(κ, λ)

    return abs(B_W(κ, λ) / 4c^2) * (2σ + 1) * 2σ * C_P(κ, λ, ξ₁) /
           abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 4)
end

function C_I_P_2_5(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    _, _, c = _abc(κ, λ)

    return abs(B_W(κ, λ) / 4c^2) * (2σ + 1) * C_P(κ, λ, ξ₁) /
           abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 4)
end

function C_I_E_dξ(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    return C_J_E(κ, ξ₁, λ)
end

function C_I_P_dξ(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    return C_J_P(κ, ξ₁, λ)
end

function C_I_E_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    return C_J_E_dκ(κ, ξ₁, λ) / ((2λ.σ + 1) * v - 2)^2
end

function C_I_E_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    return (2λ.σ + 1) * C_J_E(κ, ξ₁, λ) / abs((2λ.σ + 1) * v - 2)
end

function C_I_E_dκ_3(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    return C_J_E_dκ(κ, ξ₁, λ) / abs((2λ.σ + 1) * v - 2)
end

function C_I_P_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    _, _, c = _abc(κ, λ)

    return C_D(κ, ξ₁, λ) / abs(2c)
end

function C_I_P_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; d) = λ

    _, _, c = _abc(κ, λ)

    return C_D_dξ(κ, ξ₁, λ) + d * C_D(κ, ξ₁, λ) / abs(2c)^2
end

function C_I_P_dκ_3(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    return (
        C_D_dξ_dξ(κ, ξ₁, λ) + (2d - 1) * C_D_dξ_dξ(κ, ξ₁, λ) + d * (d - 2) * C_D(κ, ξ₁, λ)
    ) / abs((2σ + 1) * v - 2 / σ + d - 4) / abs(2c)^2
end

function C_I_P_dκ_4(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    _, _, c = _abc(κ, λ)

    return (2σ + 1) * C_D(κ, ξ₁, λ) / abs(2c)^2
end

function C_I_P_dκ_5(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < 0

    return (2σ + 1) * (2C_D_dξ(κ, ξ₁, λ) + (2d - 1) * C_D(κ, ξ₁, λ)) /
           abs((2σ + 1) * v - 2 / σ + d - 3) / abs(2c)^2
end

function C_I_P_dκ_6(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0

    return (2σ + 1) * 2σ * C_D(κ, ξ₁, λ) / abs((2σ + 1) * v - 2 / σ + d - 2) / abs(2c)^2
end

function C_I_P_dκ_7(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0

    return (2σ + 1) * C_D(κ, ξ₁, λ) / abs((2σ + 1) * v - 2 / σ + d - 2) / abs(2c)^2
end

function C_I_P_dκ_2_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0

    _, _, _, c, c_dκ = _abc_dκ(κ, λ)

    abs_BW = abs(B_W(κ, λ))
    abs_BW_dκ = abs(B_W_dκ(κ, λ))

    return (2σ + 1) * C_J_P(κ, ξ₁, λ) / abs((2σ + 1) * v - 2 / σ + d - 2)
end

function C_I_E_dξ_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    return C_J_E_dκ(κ, ξ₁, λ)
end

function C_I_E_dξ_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    return (2λ.σ + 1) * C_J_E(κ, ξ₁, λ)
end

function C_I_P_dξ_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    return C_J_P_dκ(κ, ξ₁, λ)
end

function C_I_P_dξ_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    return (2λ.σ + 1) * C_J_P(κ, ξ₁, λ)
end
