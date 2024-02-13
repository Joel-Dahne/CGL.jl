struct FunctionBounds
    P::Arb
    P_dξ::Arb
    P_dξ_dξ::Arb
    P_dξ_dξ_dξ::Arb
    P_dκ::Arb
    P_dξ_dκ::Arb
    P_dξ_dξ_dκ::Arb
    P_dϵ::Arb
    P_dξ_dϵ::Arb
    P_dξ_dξ_dϵ::Arb
    E::Arb
    E_dξ::Arb
    E_dξ_dξ::Arb
    E_dξ_dξ_dξ::Arb
    E_dκ::Arb
    E_dξ_dκ::Arb
    E_dϵ::Arb
    E_dξ_dϵ::Arb
    J_P::Arb
    J_P_dξ::Arb
    J_P_dξ_dξ::Arb
    J_P_dκ::Arb
    J_P_dϵ::Arb
    J_E::Arb
    J_E_dξ::Arb
    J_E_dξ_dξ::Arb
    J_E_dκ::Arb
    J_E_dϵ::Arb
    D::Arb
    D_dξ::Arb
    D_dξ_dξ::Arb
    H::Arb
    H_dξ::Arb
    H_dξ_dξ::Arb

    function FunctionBounds(
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb},
        CU::UBounds;
        include_dκ::Bool = false,
        include_dϵ::Bool = false,
    )
        C = new(
            C_P(κ, ϵ, ξ₁, λ, CU),
            C_P_dξ(κ, ϵ, ξ₁, λ, CU),
            C_P_dξ_dξ(κ, ϵ, ξ₁, λ),
            C_P_dξ_dξ_dξ(κ, ϵ, ξ₁, λ),
            include_dκ ? C_P_dκ(κ, ϵ, ξ₁, λ, CU) : indeterminate(κ),
            include_dκ ? C_P_dξ_dκ(κ, ϵ, ξ₁, λ, CU) : indeterminate(κ),
            include_dκ ? C_P_dξ_dξ_dκ(κ, ϵ, ξ₁, λ, CU) : indeterminate(κ),
            include_dϵ ? C_P_dϵ(κ, ϵ, ξ₁, λ, CU) : indeterminate(κ),
            include_dϵ ? C_P_dξ_dϵ(κ, ϵ, ξ₁, λ, CU) : indeterminate(κ),
            include_dϵ ? C_P_dξ_dξ_dϵ(κ, ϵ, ξ₁, λ, CU) : indeterminate(κ),
            C_E(κ, ϵ, ξ₁, λ, CU),
            C_E_dξ(κ, ϵ, ξ₁, λ, CU),
            C_E_dξ_dξ(κ, ϵ, ξ₁, λ, CU),
            C_E_dξ_dξ_dξ(κ, ϵ, ξ₁, λ, CU),
            include_dκ ? C_E_dκ(κ, ϵ, ξ₁, λ, CU) : indeterminate(κ),
            include_dκ ? C_E_dξ_dκ(κ, ϵ, ξ₁, λ, CU) : indeterminate(κ),
            include_dϵ ? C_E_dϵ(κ, ϵ, ξ₁, λ, CU) : indeterminate(κ),
            include_dϵ ? C_E_dξ_dϵ(κ, ϵ, ξ₁, λ, CU) : indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
        )

        BW = abs(B_W(κ, ϵ, λ))
        if include_dκ
            BW_dκ = abs(B_W_dκ(κ, ϵ, λ))
        end
        if include_dϵ
            BW_dϵ = abs(B_W_dϵ(κ, ϵ, λ))
        end

        C.J_P[] = C_J_P(κ, ϵ, ξ₁, λ, C, BW)
        C.J_P_dξ[] = C_J_P_dξ(κ, ϵ, ξ₁, λ, C, BW)
        C.J_P_dξ_dξ[] = C_J_P_dξ_dξ(κ, ϵ, ξ₁, λ, C, BW)
        if include_dκ
            C.J_P_dκ[] = C_J_P_dκ(κ, ϵ, ξ₁, λ, C, BW, BW_dκ)
        end
        if include_dϵ
            C.J_P_dϵ[] = C_J_P_dϵ(κ, ϵ, ξ₁, λ, C, BW, BW_dϵ)
        end

        C.J_E[] = C_J_E(κ, ϵ, ξ₁, λ, C, BW)
        C.J_E_dξ[] = C_J_E_dξ(κ, ϵ, ξ₁, λ, BW)
        C.J_E_dξ_dξ[] = C_J_E_dξ_dξ(κ, ϵ, ξ₁, λ, BW)
        if include_dκ
            C.J_E_dκ[] = C_J_E_dκ(κ, ϵ, ξ₁, λ, CU, BW, BW_dκ)
        end
        if include_dϵ
            C.J_E_dϵ[] = C_J_E_dϵ(κ, ϵ, ξ₁, λ, CU, BW, BW_dϵ)
        end

        if include_dκ
            C.D[] = C_D(κ, ϵ, ξ₁, λ, C, BW, BW_dκ)
            C.D_dξ[] = C_D_dξ(κ, ϵ, ξ₁, λ, C, BW, BW_dκ)
            C.D_dξ_dξ[] = C_D_dξ_dξ(κ, ϵ, ξ₁, λ, C, BW, BW_dκ)
        end

        if include_dϵ
            C.H[] = C_H(κ, ϵ, ξ₁, λ, C, BW, BW_dϵ)
            C.H_dξ[] = C_H_dξ(κ, ϵ, ξ₁, λ, C, BW, BW_dϵ)
            C.H_dξ_dξ[] = C_H_dξ_dξ(κ, ϵ, ξ₁, λ, C, BW, BW_dϵ)
        end

        return C
    end
end

function C_P(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return CU.U_a_b * abs(c^-a)
end

function C_E(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return CU.U_bma_b * abs((-c)^(a - b))
end

function C_P_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return abs(2c^-a) * CU.U_dz_a_b
end

function C_P_dξ_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, λ)
    z₁ = c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs(2(a + 1) * p_U(k, a + 2, b + 2, z₁) - p_U(k, a + 1, b + 1, z₁))
    end

    R = 2abs(a + 1) * C_R_U(n, a + 2, b + 2, z₁) + C_R_U(n, a + 1, b + 1, z₁)

    return abs(2a * c^-a) * (S + R * abs(z₁)^-n)
end

function C_P_dξ_dξ_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, λ)
    z₁ = c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs(-2(a + 2) * p_U(k, a + 3, b + 3, z₁) + 3p_U(k, a + 2, b + 2, z₁))
    end

    R = 2abs(a + 2) * C_R_U(n, a + 3, b + 3, z₁) + 3C_R_U(n, a + 2, b + 2, z₁)

    return abs(4a * (a + 1) * c^-a) * (S + R * abs(z₁)^-n)
end

function C_E_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    C1 = abs((-c)^(a - b)) * CU.U_bma_b
    C2 = abs((-c)^(a - b - 1)) * CU.U_dz_bma_b

    return abs(2c) * C1 + abs(2c) * C2 * ξ₁^-2
end

function C_E_dξ_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    C1 = abs((-c)^(a - b)) * CU.U_bma_b
    C2 = abs((-c)^(a - b - 1)) * CU.U_dz_bma_b
    C3 = abs((-c)^(a - b - 2)) * CU.U_dz_dz_bma_b

    return abs(4c^2 + 2c * ξ₁^-2) * C1 +
           abs(8c^2 + 2c * ξ₁^-2) * C2 * ξ₁^-2 +
           abs(4c^2) * C3 * ξ₁^-4
end

function C_E_dξ_dξ_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    C1 = abs((-c)^(a - b)) * CU.U_bma_b
    C2 = abs((-c)^(a - b - 1)) * CU.U_dz_bma_b
    C3 = abs((-c)^(a - b - 2)) * CU.U_dz_dz_bma_b
    C4 = abs((-c)^(a - b - 3)) * CU.U_dz_dz_dz_bma_b

    return abs(8c^3 + 12c^2 * ξ₁^-2) * C1 +
           abs(24c^3 + 24c^2 * ξ₁^-2) * C2 * ξ₁^-2 +
           abs(24c^3 + 12c^2 * ξ₁^-2) * C3 * ξ₁^-4 +
           abs(8c^3) * C4 * ξ₁^-6
end

# IMPROVE: This upper bound can be improved
function C_P_dκ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    C1 = CU.U_dz_a_b * abs(c^(-a - 1) * c_dκ)

    C2 = CU.U_da_a_b * abs(c^-a * a_dκ) * (2 + abs(log(c)) / log(ξ₁))

    return C1 / log(ξ₁) + C2
end

function C_E_dκ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    C1 = CU.U_bma_b * abs((-c)^(a - b) * c_dκ)

    C2 = CU.U_da_bma_b * abs((-c)^(a - b) * a_dκ) * (2 + abs(log(c)) / log(ξ₁))

    C3 = CU.U_dz_bma_b * abs((-c)^(a - b - 1) * c_dκ)

    return C1 + C2 * log(ξ₁) * ξ₁^-2 + C3 * ξ₁^-2
end

# IMPROVE: This upper bound can be improved
function C_P_dξ_dκ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    C1 = CU.U_dz_a_b * abs(c^(-a - 1) * 2c_dκ)

    C2 = CU.U_dz_dz_a_b * abs(c^(-a - 1) * 2c_dκ)

    C3 = CU.U_da_ap1_bp1 * abs(c^-a * 2a_dκ * a) * (2 + abs(log(c)) / log(ξ₁))

    C4 = CU.U_ap1_bp1 * abs(c^-a * 2a_dκ)

    return C1 / log(ξ₁) + C2 / log(ξ₁) + C3 + C4 / log(ξ₁)
end

function C_E_dξ_dκ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    C1 = 2CU.U_bma_b * abs((-c)^(a - b) * c_dκ) * (abs(c) + ξ₁^-2)

    C2 = 2CU.U_dz_bma_b * abs((-c)^(a - b - 1) * c_dκ) * (abs(2c) + ξ₁^-2)

    C3 = 2CU.U_da_bma_b * abs((-c)^(a - b) * a_dκ * c) * (2 + abs(log(-c)) / log(ξ₁))

    C4 = 2CU.U_dz_dz_bma_b * abs((-c)^(a - b - 2) * c * c_dκ)

    C5 = 2CU.U_bmap1_bp1 * abs((-c)^(a - b - 1) * a_dκ * c)

    C6 =
        2CU.U_da_bmap1_bp1 *
        abs((b - a) * (-c)^(a - b - 1) * a_dκ * c) *
        (2 + abs(log(-c)) / log(ξ₁))

    return C1 +
           C2 * ξ₁^-2 +
           C3 * log(ξ₁) * ξ₁^-2 +
           C4 * ξ₁^-2 +
           C5 * ξ₁^-4 +
           C6 * log(ξ₁) * ξ₁^-4
end

# IMPROVE: This upper bound can be improved
function C_P_dξ_dξ_dκ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    C1 = CU.U_dz_a_b * abs(c^(-a - 1) * 2c_dκ)

    C2 = CU.U_dz_dz_a_b * abs(c^(-a - 1) * 10c_dκ)

    C3 = CU.U_dz_dz_dz_a_b * abs(c^(-a - 1) * 4c_dκ)

    C4 = CU.U_ap1_bp1 * abs(c^-a * 2a_dκ)

    C5 = CU.U_da_ap1_bp1 * abs(c^-a * 2a_dκ * a) * (2 + abs(log(c)) / log(ξ₁))

    C6 = CU.U_ap2_bp2 * abs(c^-a * 4a_dκ * (2a + 1))

    C7 = CU.U_da_ap2_bp2 * abs(c^-a * 4a_dκ * a * (a + 1)) * (2 + abs(log(c)) / log(ξ₁))

    return C1 / log(ξ₁) +
           C2 / log(ξ₁) +
           C3 / log(ξ₁) +
           C4 / log(ξ₁) +
           C5 +
           C6 / log(ξ₁) +
           C7
end

function C_P_dϵ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    C = CU.U_dz_a_b * abs(c^(-a - 1) * c_dϵ)

    return C
end

function C_E_dϵ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    C1 = CU.U_bma_b * abs((-c)^(a - b) * c_dϵ)

    C2 = CU.U_dz_bma_b * abs((-c)^(a - b - 1) * c_dϵ)

    return C1 + C2 * ξ₁^-2
end

# IMPROVE: This upper bound can be improved
function C_P_dξ_dϵ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    C1 = CU.U_dz_a_b * abs(c^(-a - 1) * 2c_dϵ)

    C2 = CU.U_dz_dz_a_b * abs(c^(-a - 1) * 2c_dϵ)

    return C1 + C2
end

function C_E_dξ_dϵ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    C1 = 2CU.U_bma_b * abs((-c)^(a - b) * c_dϵ) * (abs(c) + ξ₁^-2)

    C2 = 2CU.U_dz_bma_b * abs((-c)^(a - b - 1) * c_dϵ) * (abs(2c) + ξ₁^-2)

    C3 = 2CU.U_dz_dz_bma_b * abs((-c)^(a - b - 2) * c * c_dϵ)

    return C1 + C2 * ξ₁^-2 + C3 * ξ₁^-2
end

# IMPROVE: This upper bound can be improved
function C_P_dξ_dξ_dϵ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    C1 = CU.U_dz_a_b * abs(c^(-a - 1) * 2c_dϵ)

    C2 = CU.U_dz_dz_a_b * abs(c^(-a - 1) * 10c_dϵ)

    C3 = CU.U_dz_dz_dz_a_b * abs(c^(-a - 1) * 4c_dϵ)

    return C1 + C2 + C3
end

C_J_P(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C::FunctionBounds, BW::Arb) = BW * C.P

C_J_E(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C::FunctionBounds, BW::Arb) = BW * C.E

function C_J_P_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C::FunctionBounds, BW::Arb)
    _, _, c = _abc(κ, ϵ, λ)
    (; d) = λ

    return BW * (C.P * (abs(2c) + (d - 1) * ξ₁^-2) + C.P_dξ * ξ₁^-2)
end

function C_J_E_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, BW::Arb)
    a, b, c = _abc(κ, ϵ, λ)
    (; d) = λ
    z₁ = -c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs((d - 1) * p_U(k, b - a, b, z₁) - 2(b - a) * p_U(k, b - a + 1, b + 1, z₁))
    end

    R = (d - 1) * C_R_U(n, b - a, b, z₁) + abs(2(b - a)) * C_R_U(n, b - a + 1, b + 1, z₁)

    return BW * abs((-c)^(a - b)) * (S + R * abs(z₁)^-n)
end

function C_J_P_dξ_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C::FunctionBounds, BW::Arb)
    _, _, c = _abc(κ, ϵ, λ)
    (; d) = λ

    return BW * (
        C.P * (abs(4c^2) + abs(2c) * (2d - 1) * ξ₁^-2 + (d - 1) * (d - 2) * ξ₁^-4) +
        C.P_dξ * (abs(4c) + 2(d - 1) * ξ₁^-2) * ξ₁^-2 +
        C.P_dξ_dξ * ξ₁^-4
    )
end

function C_J_E_dξ_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, BW::Arb)
    a, b, c = _abc(κ, ϵ, λ)
    (; d) = λ
    z₁ = -c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs(
            (d - 1) * (d - 2) * p_U(k, b - a, b, z₁) -
            2(2d - 1) * (b - a) * p_U(k, b - a + 1, b + 1, z₁) +
            4(b - a) * (b - a + 1) * p_U(k, b - a + 2, b + 2, z₁),
        )
    end

    R =
        (d - 1) * (d - 2) * C_R_U(n, b - a, b, z₁) +
        abs(2(2d - 1) * (b - a)) * C_R_U(n, b - a + 1, b + 1, z₁) +
        abs(4(b - a) * (b - a + 1)) * C_R_U(n, b - a + 2, b + 2, z₁)

    return BW * abs((-c)^(a - b)) * (S + R * abs(z₁)^-n)
end

function C_J_P_dκ(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    BW::Arb,
    BW_dκ::Arb,
)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    C1 = C.P * (abs(c_dκ) * BW + BW_dκ * ξ₁^-2)

    C2 = C.P_dκ * BW

    return C1 + C2 * log(ξ₁) * ξ₁^-2
end

# IMPROVE: This upper bound can be improved
function C_J_E_dκ(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    CU::UBounds,
    BW::Arb,
    BW_dκ::Arb,
)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    C1 = BW * abs((-c)^(a - b) * a_dκ) * (2 + abs(log(-c)) / log(ξ₁)) * CU.U_da_bma_b

    C2 = BW_dκ * abs((-c)^(a - b)) * CU.U_bma_b

    C3 = BW * abs((-c)^(a - b - 1) * c_dκ) * CU.U_dz_bma_b

    return C1 + C2 / log(ξ₁) + C3 / log(ξ₁)
end

function C_J_P_dϵ(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    BW::Arb,
    BW_dϵ::Arb,
)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    C1 = C.P * (abs(c_dϵ) * BW + BW_dϵ * ξ₁^-2)

    C2 = C.P_dϵ * BW

    return C1 + C2 * ξ₁^-2
end

# IMPROVE: This upper bound can be improved
function C_J_E_dϵ(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    CU::UBounds,
    BW::Arb,
    BW_dϵ::Arb,
)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    C1 = BW_dϵ * abs((-c)^(a - b)) * CU.U_bma_b

    C2 = BW * abs((-c)^(a - b - 1) * c_dϵ) * CU.U_dz_bma_b

    return C1 + C2
end

function C_D(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    BW::Arb,
    BW_dκ::Arb,
)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    C1 = abs(c_dκ) * BW * C.P

    C2 = BW_dκ * C.P

    C3 = BW * C.P_dκ

    return C1 + (C2 + C3 * log(ξ₁)) * ξ₁^-2
end

function C_D_dξ(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    BW::Arb,
    BW_dκ::Arb,
)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    C1 = abs(c_dκ) * BW * C.P_dξ

    C2 = BW_dκ * C.P_dξ

    C3 = 2BW_dκ * C.P

    C4 = BW * C.P_dξ_dκ

    C5 = 2BW * C.P_dκ

    return C1 + (C2 + C3 + (C4 + C5) * log(ξ₁)) * ξ₁^-2
end

function C_D_dξ_dξ(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    BW::Arb,
    BW_dκ::Arb,
)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    C1 = abs(c_dκ) * BW * C.P_dξ_dξ

    C2 = BW_dκ * C.P_dξ_dξ

    C3 = 4BW_dκ * C.P_dξ

    C4 = 6BW_dκ * C.P

    C5 = BW * C.P_dξ_dξ_dκ

    C6 = 2BW * C.P_dξ_dκ

    C7 = 6BW * C.P_dκ

    return C1 + (C2 + C3 + C4 + (C5 + C6 + C7) * log(ξ₁)) * ξ₁^-2
end

function C_H(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    BW::Arb,
    BW_dϵ::Arb,
)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    C1 = abs(c_dϵ) * BW * C.P

    C2 = BW_dϵ * C.P

    C3 = BW * C.P_dϵ

    return C1 + (C2 + C3) * ξ₁^-2
end

function C_H_dξ(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    BW::Arb,
    BW_dϵ::Arb,
)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    C1 = abs(c_dϵ) * BW * C.P_dξ

    C2 = BW_dϵ * C.P_dξ

    C3 = 2BW_dϵ * C.P

    C4 = BW * C.P_dξ_dϵ

    C5 = 2BW * C.P_dϵ

    return C1 + (C2 + C3 + C4 + C5) * ξ₁^-2
end

function C_H_dξ_dξ(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    BW::Arb,
    BW_dϵ::Arb,
)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    C1 = abs(c_dϵ) * BW * C.P_dξ_dξ

    C2 = BW_dϵ * C.P_dξ_dξ

    C3 = 4BW_dϵ * C.P_dξ

    C4 = 6BW_dϵ * C.P

    C5 = BW * C.P_dξ_dξ_dϵ

    C6 = 2BW * C.P_dξ_dϵ

    C7 = 6BW * C.P_dϵ

    return C1 + (C2 + C3 + C4 + C5 + C6 + C7) * ξ₁^-2
end
