function C_I_E(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_J_E(κ, ξ₁, λ) / abs((2λ.σ + 1) * v - 2)
end

function C_I_P(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_J_P(κ, ξ₁, λ) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 2)
end

function C_I_P_1(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    (; σ, d) = λ

    _, _, c = _abc(κ, λ)

    return C_J_P(κ, ξ₁, λ) / abs(2c) +
           abs(B_W(κ, λ) / 2c) * (
        C_P_dξ(κ, λ, ξ₁) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 3) +
        abs(d - 2) * C_P(κ, λ, ξ₁) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 4)
    )
end

function C_I_P_2(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    (; σ, d) = λ

    _, _, c = _abc(κ, λ)

    return abs(B_W(κ, λ) / 2c) * (2σ + 1) * C_P(κ, λ, ξ₁) /
           abs((2σ + 1) * v - 2 / σ + d - 3)
end

function C_I_E_dξ(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_J_E(κ, ξ₁, λ)
end

function C_I_P_dξ(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_J_P(κ, ξ₁, λ)
end

function C_I_E_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_J_E_dκ(κ, ξ₁, λ) / ((2λ.σ + 1) * v - 2)^2
end

function C_I_E_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return (2λ.σ + 1) * C_J_E(κ, ξ₁, λ) / abs((2λ.σ + 1) * v - 2)
end

function C_I_E_dκ_3(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_J_E_dκ(κ, ξ₁, λ) / abs((2λ.σ + 1) * v - 2)
end

function C_I_P_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    _, _, _, c, c_dκ = _abc_dκ(κ, λ)

    abs_BW = abs(B_W(κ, λ))
    abs_BW_dκ = abs(B_W_dκ(κ, λ))

    return inv(abs(2c)) * (
        abs(c_dκ * abs_BW) * C_P(κ, λ, ξ₁) +
        abs_BW * C_P_dκ(κ, λ, ξ₁) * log(ξ₁) * ξ₁^(-2) +
        abs_BW_dκ * C_P(κ, λ, ξ₁) * ξ₁^(-2)
    )
end

function C_I_P_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    (; σ, d) = λ

    _, _, _, c, c_dκ = _abc_dκ(κ, λ)

    abs_BW = abs(B_W(κ, λ))
    abs_BW_dκ = abs(B_W_dκ(κ, λ))

    return inv(abs(2c * ((2σ + 1) * v - 2 / σ + d - 2))) * (
        abs(c_dκ * abs_BW) * (C_P_dξ(κ, λ, ξ₁) + C_P(κ, λ, ξ₁) * (2 + abs(d - 2))) +
        abs_BW * (C_P_dξ_dκ(κ, λ, ξ₁) + C_P_dκ(κ, λ, ξ₁) * abs(d - 2)) * log(ξ₁) * ξ₁^(-2) +
        abs_BW_dκ * (C_P_dξ(κ, λ, ξ₁) + C_P(κ, λ, ξ₁) * abs(d - 2)) * ξ₁^(-2)
    )
end

function C_I_P_dκ_3(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    (; σ, d) = λ

    _, _, _, c, c_dκ = _abc_dκ(κ, λ)

    abs_BW = abs(B_W(κ, λ))
    abs_BW_dκ = abs(B_W_dκ(κ, λ))

    return (2σ + 1) * (
        abs(c_dκ * abs_BW) * C_P(κ, λ, ξ₁) +
        abs_BW * C_P_dκ(κ, λ, ξ₁) * log(ξ₁) * ξ₁^(-2) +
        abs_BW_dκ * C_P(κ, λ, ξ₁) * ξ₁^(-2)
    ) / abs(2c * ((2σ + 1) * v - 2 / σ + d - 1))
end

function C_I_P_dκ_4(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    (; σ, d) = λ

    _, _, _, c, c_dκ = _abc_dκ(κ, λ)

    abs_BW = abs(B_W(κ, λ))
    abs_BW_dκ = abs(B_W_dκ(κ, λ))

    return (2σ + 1) * C_J_P(κ, ξ₁, λ) / abs((2σ + 1) * v - 2 / σ + d - 2)
end

function C_I_E_dξ_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_J_E_dκ(κ, ξ₁, λ)
end

function C_I_E_dξ_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return (2λ.σ + 1) * C_J_E(κ, ξ₁, λ)
end

function C_I_P_dξ_dκ_1(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return C_J_P_dκ(κ, ξ₁, λ)
end

function C_I_P_dξ_dκ_2(κ::Arb, ξ₁::Arb, v::Arb, λ::AbstractGLParams{Arb})
    return (2λ.σ + 1) * C_J_P(κ, ξ₁, λ)
end
