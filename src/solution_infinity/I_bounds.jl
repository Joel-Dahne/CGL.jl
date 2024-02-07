function C_I_E(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    @assert (2λ.σ + 1) * v - 2 < 0
    return C.J_E / abs((2λ.σ + 1) * v - 2)
end

function C_I_P(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    @assert (2λ.σ + 1) * v - 2 / λ.σ + λ.d - 2 < 0
    return C.J_P / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 2)
end

function C_I_P_1_1(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 3 < 0

    _, _, c = _abc(κ, ϵ, λ)

    return abs(B_W(κ, ϵ, λ) / 2c) * (
        C.P +
        C.P_dξ / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 4) +
        abs(d - 2) * C.P / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 4)
    )
end

function C_I_P_1_2(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    @assert (2λ.σ + 1) * v - 2 / λ.σ + λ.d - 3 < 0

    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    return abs(B_W(κ, ϵ, λ) / 2c) * (2σ + 1) * C.P / abs((2σ + 1) * v - 2 / σ + d - 3)
end

function C_I_P_2_1(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    _, _, c = _abc(κ, ϵ, λ)

    return abs(B_W(κ, ϵ, λ) / 2c) * C.P
end

function C_I_P_2_2(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    _, _, c = _abc(κ, ϵ, λ)

    return abs(B_W(κ, ϵ, λ) / 4c^2) * (
        C.P_dξ +
        abs(d - 2) * C.P +
        (C.P_dξ_dξ + abs(2d - 5) * C.P_dξ + abs((d - 2) * (d - 4)) * C.P) /
        abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 6)
    )
end

function C_I_P_2_3(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    _, _, c = _abc(κ, ϵ, λ)

    return abs(B_W(κ, ϵ, λ) / 4c^2) *
           (2σ + 1) *
           (C.P + (2C.P_dξ + abs(2d - 5) * C.P) / abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 5))
end

function C_I_P_2_4(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    _, _, c = _abc(κ, ϵ, λ)

    return abs(B_W(κ, ϵ, λ) / 4c^2) * (2σ + 1) * 2σ * C.P /
           abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 4)
end

function C_I_P_2_5(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    _, _, c = _abc(κ, ϵ, λ)

    return abs(B_W(κ, ϵ, λ) / 4c^2) * (2σ + 1) * C.P /
           abs((2λ.σ + 1) * v - 2 / λ.σ + λ.d - 4)
end

function C_I_E_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.J_E
end

function C_I_P_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.J_P
end

function C_I_E_dκ(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.J_E_dκ * (inv(((2λ.σ + 1) * v - 2)^2) + log(ξ₁) / abs((2λ.σ + 1) * v - 2))
end

function C_I_P_dκ_1_1(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    _, _, c = _abc(κ, ϵ, λ)

    return C.D / abs(2c)
end

function C_I_P_dκ_1_2(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    return C.D_dξ + d * C.D / abs(2c)^2
end

function C_I_P_dκ_1_3(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    return (C.D_dξ_dξ + abs(2d - 1) * C.D_dξ + abs(d * (d - 2)) * C.D) /
           abs((2σ + 1) * v - 2 / σ + d - 4) / abs(2c)^2
end

function C_I_P_dκ_1_4(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    return (2σ + 1) * C.D / abs(2c)^2
end

function C_I_P_dκ_1_5(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < 0

    return (2σ + 1) * (2C.D_dξ + abs(2d - 1) * C.D) / abs((2σ + 1) * v - 2 / σ + d - 3) /
           abs(2c)^2
end

function C_I_P_dκ_1_6(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0

    return (2σ + 1) * 2σ * C.D / abs((2σ + 1) * v - 2 / σ + d - 2) / abs(2c)^2
end

function C_I_P_dκ_1_7(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0

    return (2σ + 1) * C.D / abs((2σ + 1) * v - 2 / σ + d - 2) / abs(2c)^2
end

function C_I_E_dϵ(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    return C.J_E_dϵ / abs((2λ.σ + 1) * v - 2)
end

function C_I_P_dϵ_1_1(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    _, _, c = _abc(κ, ϵ, λ)

    return C.H / abs(2c)
end

function C_I_P_dϵ_1_2(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    return C.H_dξ + d * C.H / abs(2c)^2
end

function C_I_P_dϵ_1_3(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0

    return (C.H_dξ_dξ + abs(2d - 1) * C.H_dξ + abs(d * (d - 2)) * C.H) /
           abs((2σ + 1) * v - 2 / σ + d - 4) / abs(2c)^2
end

function C_I_P_dϵ_1_4(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    return (2σ + 1) * C.H / abs(2c)^2
end

function C_I_P_dϵ_1_5(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < 0

    return (2σ + 1) * (2C.H_dξ + abs(2d - 1) * C.H) / abs((2σ + 1) * v - 2 / σ + d - 3) /
           abs(2c)^2
end

function C_I_P_dϵ_1_6(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0

    return (2σ + 1) * 2σ * C.H / abs((2σ + 1) * v - 2 / σ + d - 2) / abs(2c)^2
end

function C_I_P_dϵ_1_7(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = λ

    _, _, c = _abc(κ, ϵ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0

    return (2σ + 1) * C.H / abs((2σ + 1) * v - 2 / σ + d - 2) / abs(2c)^2
end
