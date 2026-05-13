"""
    FunctionEnclosures(κ, ϵ, ξ₁, Λ; include_dκ = false, include_dϵ = false)

Contains enclosures of the functions

- [`P](@ref)
- [`P_dξ`](@ref)
- [`P_dξ_dξ`](@ref)
- [`E`](@ref)
- [`E_dξ`](@ref)
- [`J_P`](@ref)
- [`J_E`](@ref)

when evaluated at `ξ₁`.

If `include_dκ = true` the also include enclosures for:

- [`P_dκ`](@ref)
- [`P_dξ_dκ`](@ref)
- [`E_dκ`](@ref)
- [`E_dξ_dκ`](@ref)
- [`D`](@ref)
- [`D_dξ`](@ref)
- [`J_P_dκ`](@ref)
- [`J_E_dκ`](@ref)

If `include_dϵ = true` the also include enclosures for:

- [`P_dϵ`](@ref)
- [`P_dξ_dϵ`](@ref)
- [`E_dϵ`](@ref)
- [`E_dξ_dϵ`](@ref)
- [`H`](@ref)
- [`H_dξ`](@ref)
- [`J_P_dϵ`](@ref)
- [`J_E_dϵ`](@ref)
"""
struct FunctionEnclosures
    # Always included
    P::Acb
    P_dξ::Acb
    P_dξ_dξ::Acb
    E::Acb
    E_dξ::Acb
    J_P::Acb
    J_E::Acb
    # Included when include_dκ = true (otherwise indeterminate)
    P_dκ::Acb
    P_dξ_dκ::Acb
    E_dκ::Acb
    E_dξ_dκ::Acb
    D::Acb
    D_dξ::Acb
    J_P_dκ::Acb
    J_E_dκ::Acb
    # Included when include_dϵ = true (otherwise indeterminate)
    P_dϵ::Acb
    P_dξ_dϵ::Acb
    E_dϵ::Acb
    E_dξ_dϵ::Acb
    H::Acb
    H_dξ::Acb
    J_P_dϵ::Acb
    J_E_dϵ::Acb

    FunctionEnclosures() = new(
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
    )
end

function FunctionEnclosures(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb};
    include_dκ::Bool = false,
    include_dϵ::Bool = false,
)
    F = FunctionEnclosures()

    ξ₁_series = ArbSeries((ξ₁, 1))
    BW = B_W(κ, ϵ, Λ)

    P_series = P(ArbSeries(ξ₁_series, degree = 2), κ, ϵ, Λ)
    F.P[] = P_series[0]
    F.P_dξ[] = P_series[1]
    F.P_dξ_dξ[] = 2P_series[2]

    E_series = E(ξ₁_series, κ, ϵ, Λ)
    F.E[] = E_series[0]
    F.E_dξ[] = E_series[1]

    F.J_P[] = J_P(ξ₁, κ, ϵ, Λ, p = F.P)
    F.J_E[] = J_E(ξ₁, κ, ϵ, Λ, e = F.E)

    if include_dκ
        BW_dκ = B_W_dκ(κ, ϵ, Λ)

        P_dκ_series = P_dκ(ξ₁_series, κ, ϵ, Λ)
        F.P_dκ[] = P_dκ_series[0]
        F.P_dξ_dκ[] = P_dκ_series[1]

        E_dκ_series = E_dκ(ξ₁_series, κ, ϵ, Λ)
        F.E_dκ[] = E_dκ_series[0]
        F.E_dξ_dκ[] = E_dκ_series[1]

        D_series = D(ξ₁_series, κ, ϵ, Λ; p = P_series, p_dκ = P_dκ_series, BW, BW_dκ)
        F.D[] = D_series[0]
        F.D_dξ[] = D_series[1]

        F.J_P_dκ[] = J_P_dκ(ξ₁, κ, ϵ, Λ, d = F.D)
        F.J_E_dκ[] = J_E_dκ(ξ₁, κ, ϵ, Λ, e = F.E, e_dκ = F.E_dκ; BW, BW_dκ)
    end

    if include_dϵ
        BW_dϵ = B_W_dϵ(κ, ϵ, Λ)

        P_dϵ_series = P_dϵ(ξ₁_series, κ, ϵ, Λ)
        F.P_dϵ[] = P_dϵ_series[0]
        F.P_dξ_dϵ[] = P_dϵ_series[1]

        E_dϵ_series = E_dϵ(ξ₁_series, κ, ϵ, Λ)
        F.E_dϵ[] = E_dϵ_series[0]
        F.E_dξ_dϵ[] = E_dϵ_series[1]

        H_series = H(ξ₁_series, κ, ϵ, Λ; p = P_series, p_dϵ = P_dϵ_series, BW, BW_dϵ)
        F.H[] = H_series[0]
        F.H_dξ[] = H_series[1]

        F.J_P_dϵ[] = J_P_dϵ(ξ₁, κ, ϵ, Λ, h = F.H)
        F.J_E_dϵ[] = J_E_dϵ(ξ₁, κ, ϵ, Λ; e = F.E, e_dϵ = F.E_dϵ, BW, BW_dϵ)
    end

    return F
end

function P(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)

    z = c * ξ^2

    return U(a, b, z)
end

function P_dξ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ

    return U_dz(a, b, z) * z_dξ
end

function P_dξ_dξ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c

    return U_dz(a, b, z, 2) * z_dξ^2 + U_dz(a, b, z) * z_dξ_dξ
end

function P_dξ_dξ_dξ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c
    # z_dξ_dξ_dξ = 0

    return U_dz(a, b, z, 3) * z_dξ^3 + U_dz(a, b, z, 2) * 3z_dξ * z_dξ_dξ
end

function P_dκ(ξ, κ, ϵ, Λ::CGLParams)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, Λ)

    z = c * ξ^2
    z_dκ = c_dκ * ξ^2

    return U_da(a, b, z) * a_dκ + U_dz(a, b, z) * z_dκ
end

function P_dξ_dκ(ξ, κ, ϵ, Λ::CGLParams)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dκ = c_dκ * ξ^2
    z_dξ_dκ = 2c_dκ * ξ

    return (U_dzda(a, b, z) * a_dκ + U_dz(a, b, z, 2) * z_dκ) * z_dξ +
           U_dz(a, b, z) * z_dξ_dκ
end

function P_dξ_dξ_dκ(ξ, κ, ϵ, Λ::CGLParams)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c
    z_dκ = c_dκ * ξ^2
    z_dξ_dκ = 2c_dκ * ξ
    z_dξ_dξ_dκ = 2c_dκ

    return (U_dzda(a, b, z, 2) * a_dκ + U_dz(a, b, z, 3) * z_dκ) * z_dξ^2 +
           U_dz(a, b, z, 2) * 2z_dξ * z_dξ_dκ +
           (U_dzda(a, b, z) * a_dκ + U_dz(a, b, z, 2) * z_dκ) * z_dξ_dξ +
           U_dz(a, b, z) * z_dξ_dξ_dκ
end

function P_dϵ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, Λ)

    z = c * ξ^2
    z_dϵ = c_dϵ * ξ^2

    return U_dz(a, b, z) * z_dϵ
end

function P_dξ_dϵ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dϵ = c_dϵ * ξ^2
    z_dξ_dϵ = 2c_dϵ * ξ

    return U_dz(a, b, z, 2) * z_dϵ * z_dξ + U_dz(a, b, z) * z_dξ_dϵ
end

function P_dξ_dξ_dϵ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c
    z_dϵ = c_dϵ * ξ^2
    z_dξ_dϵ = 2c_dϵ * ξ
    z_dξ_dξ_dϵ = 2c_dϵ

    return U_dz(a, b, z, 3) * z_dϵ * z_dξ^2 +
           U_dz(a, b, z, 2) * 2z_dξ * z_dξ_dϵ +
           U_dz(a, b, z, 2) * z_dϵ * z_dξ_dξ +
           U_dz(a, b, z) * z_dξ_dξ_dϵ
end

function E(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)

    z = c * ξ^2

    return exp(z) * U(b - a, b, -z)
end

function E_dξ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ

    return exp(z) * (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ
end

function E_dξ_dξ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c

    return exp(z) * (
        U(b - a, b, -z) * (z_dξ^2 + z_dξ_dξ) - U_dz(b - a, b, -z) * (2z_dξ^2 + z_dξ_dξ) +
        U_dz(b - a, b, -z, 2) * z_dξ^2
    )
end

function E_dξ_dξ_dξ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c
    #z_dξ_dξ_dξ = 0

    return exp(z) * (
        U(b - a, b, -z) * (z_dξ^3 + 3z_dξ * z_dξ_dξ) -
        U_dz(b - a, b, -z) * (3z_dξ^3 + 6z_dξ * z_dξ_dξ) +
        U_dz(b - a, b, -z, 2) * (3z_dξ^3 + 3z_dξ * z_dξ_dξ) -
        U_dz(b - a, b, -z, 3) * z_dξ^3
    )
end

function E_dκ(ξ, κ, ϵ, Λ::CGLParams)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, Λ)

    z = c * ξ^2
    z_dκ = c_dκ * ξ^2

    return exp(z) * (
        z_dκ * U(b - a, b, -z) +
        (U_da(b - a, b, -z) * (-a_dκ) + U_dz(b - a, b, -z) * (-z_dκ))
    )
end

function E_dξ_dκ(ξ, κ, ϵ, Λ::CGLParams)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dκ = c_dκ * ξ^2
    z_dξ_dκ = 2c_dκ * ξ

    return exp(z) * (
        (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ * z_dκ +
        (
            U_da(b - a, b, -z) * (-a_dκ) + U_dz(b - a, b, -z) * (-z_dκ) -
            (U_dzda(b - a, b, -z) * (-a_dκ) + U_dz(b - a, b, -z, 2) * (-z_dκ))
        ) * z_dξ +
        (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ_dκ
    )
end

function E_dϵ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, Λ)

    z = c * ξ^2
    z_dϵ = c_dϵ * ξ^2

    return exp(z) * (z_dϵ * U(b - a, b, -z) + U_dz(b - a, b, -z) * (-z_dϵ))
end

function E_dξ_dϵ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, Λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dϵ = c_dϵ * ξ^2
    z_dξ_dϵ = 2c_dϵ * ξ

    return exp(z) * (
        (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ * z_dϵ +
        (U_dz(b - a, b, -z) * (-z_dϵ) - U_dz(b - a, b, -z, 2) * (-z_dϵ)) * z_dξ +
        (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ_dϵ
    )
end

function W(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)

    z = c * ξ^2

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 2c * exp(sgn * im * (b - a) * π) * ξ * z^-b * exp(z)
end

function J_P(ξ, κ, ϵ, Λ::CGLParams; p = P(ξ, κ, ϵ, Λ))
    c = _c(κ, ϵ, Λ)

    return B_W(κ, ϵ, Λ) * p * exp(-c * ξ^2) * ξ^(Λ.d - 1)
end

function J_E(ξ, κ, ϵ, Λ::CGLParams; e = E(ξ, κ, ϵ, Λ))
    c = _c(κ, ϵ, Λ)

    return B_W(κ, ϵ, Λ) * e * exp(-c * ξ^2) * ξ^(Λ.d - 1)
end

# These four are only used for testing and are not performance critical
J_P_dξ(ξ, κ, ϵ, Λ::CGLParams) = J_P(ArbSeries((ξ, 1)), κ, ϵ, Λ)[1]
J_E_dξ(ξ, κ, ϵ, Λ::CGLParams) = J_E(ArbSeries((ξ, 1)), κ, ϵ, Λ)[1]
J_P_dξ_dξ(ξ, κ, ϵ, Λ::CGLParams) = 2J_P(ArbSeries((ξ, 1), degree = 2), κ, ϵ, Λ)[2]
J_E_dξ_dξ(ξ, κ, ϵ, Λ::CGLParams) = 2J_E(ArbSeries((ξ, 1), degree = 2), κ, ϵ, Λ)[2]

function J_P_dκ(ξ, κ, ϵ, Λ::CGLParams; d = D(ξ, κ, ϵ, Λ))
    c = _c(κ, ϵ, Λ)

    return d * exp(-c * ξ^2) * ξ^(Λ.d + 1)
end

function J_E_dκ(
    ξ,
    κ,
    ϵ,
    Λ::CGLParams;
    e = E(ξ, κ, ϵ, Λ),
    e_dκ = E_dκ(ξ, κ, ϵ, Λ),
    BW = B_W(κ, ϵ, Λ),
    BW_dκ = B_W_dκ(κ, ϵ, Λ),
)
    c, c_dκ = _c(κ, ϵ, Λ), _c_dκ(κ, ϵ, Λ)

    return (BW_dκ * e * ξ^-2 + BW * e_dκ * ξ^-2 - BW * e * c_dκ) *
           exp(-c * ξ^2) *
           ξ^(Λ.d + 1)
end

function J_P_dϵ(ξ, κ, ϵ, Λ::CGLParams; h = H(ξ, κ, ϵ, Λ))
    c = _c(κ, ϵ, Λ)

    return h * exp(-c * ξ^2) * ξ^(Λ.d + 1)
end

function J_E_dϵ(
    ξ,
    κ,
    ϵ,
    Λ::CGLParams;
    e = E(ξ, κ, ϵ, Λ),
    e_dϵ = E_dϵ(ξ, κ, ϵ, Λ),
    BW = B_W(κ, ϵ, Λ),
    BW_dϵ = B_W_dϵ(κ, ϵ, Λ),
)
    c, c_dϵ = _c(κ, ϵ, Λ), _c_dϵ(κ, ϵ, Λ)

    return (BW_dϵ * e * ξ^-2 + BW * e_dϵ * ξ^-2 - BW * e * c_dϵ) *
           exp(-c * ξ^2) *
           ξ^(Λ.d + 1)
end

J_P_dξ(ξ::Float64, κ::Float64, ϵ::Float64, Λ::CGLParams{Float64}) =
    ComplexF64(J_P(ArbSeries((ξ, 1)), κ, ϵ, Λ)[1])
J_E_dξ(ξ::Float64, κ::Float64, ϵ::Float64, Λ::CGLParams{Float64}) =
    ComplexF64(J_E(ArbSeries((ξ, 1)), κ, ϵ, Λ)[1])
J_P_dξ_dξ(ξ::Float64, κ::Float64, ϵ::Float64, Λ::CGLParams{Float64}) =
    ComplexF64(2J_P(ArbSeries((ξ, 1), degree = 2), κ, ϵ, Λ)[2])
J_E_dξ_dξ(ξ::Float64, κ::Float64, ϵ::Float64, Λ::CGLParams{Float64}) =
    ComplexF64(2J_E(ArbSeries((ξ, 1), degree = 2), κ, ϵ, Λ)[2])
J_P_dκ(ξ::Float64, κ::Float64, ϵ::Float64, Λ::CGLParams{Float64}) =
    ComplexF64(J_P(ξ, ArbSeries((κ, 1)), ϵ, Λ)[1])
J_E_dκ(ξ::Float64, κ::Float64, ϵ::Float64, Λ::CGLParams{Float64}) =
    ComplexF64(J_E(ξ, ArbSeries((κ, 1)), ϵ, Λ)[1])
J_P_dϵ(ξ::Float64, κ::Float64, ϵ::Float64, Λ::CGLParams{Float64}) =
    ComplexF64(J_P(ξ, κ, ArbSeries((ϵ, 1)), Λ)[1])
J_E_dϵ(ξ::Float64, κ::Float64, ϵ::Float64, Λ::CGLParams{Float64}) =
    ComplexF64(J_E(ξ, κ, ArbSeries((ϵ, 1)), Λ)[1])

function D(
    ξ,
    κ,
    ϵ,
    Λ::CGLParams;
    p = P(ξ, κ, ϵ, Λ),
    p_dκ = P_dκ(ξ, κ, ϵ, Λ),
    BW = B_W(κ, ϵ, Λ),
    BW_dκ = B_W_dκ(κ, ϵ, Λ),
)
    c_dκ = _c_dκ(κ, ϵ, Λ)

    return p * (-c_dκ * BW + BW_dκ * ξ^-2) + BW * p_dκ * ξ^-2
end

function D_dξ(ξ, κ, ϵ, Λ::CGLParams)
    c_dκ = _c_dκ(κ, ϵ, Λ)

    return -c_dκ * B_W(κ, ϵ, Λ) * P_dξ(ξ, κ, ϵ, Λ) +
           B_W_dκ(κ, ϵ, Λ) * P_dξ(ξ, κ, ϵ, Λ) * ξ^-2 +
           -2B_W_dκ(κ, ϵ, Λ) * P(ξ, κ, ϵ, Λ) * ξ^-3 +
           B_W(κ, ϵ, Λ) * P_dξ_dκ(ξ, κ, ϵ, Λ) * ξ^-2 +
           -2B_W(κ, ϵ, Λ) * P_dκ(ξ, κ, ϵ, Λ) * ξ^-3
end

function D_dξ_dξ(ξ, κ, ϵ, Λ::CGLParams)
    c_dκ = _c_dκ(κ, ϵ, Λ)

    return -c_dκ * B_W(κ, ϵ, Λ) * P_dξ_dξ(ξ, κ, ϵ, Λ) +
           B_W_dκ(κ, ϵ, Λ) * P_dξ_dξ(ξ, κ, ϵ, Λ) * ξ^-2 +
           -4B_W_dκ(κ, ϵ, Λ) * P_dξ(ξ, κ, ϵ, Λ) * ξ^-3 +
           6B_W_dκ(κ, ϵ, Λ) * P(ξ, κ, ϵ, Λ) * ξ^-4 +
           B_W(κ, ϵ, Λ) * P_dξ_dξ_dκ(ξ, κ, ϵ, Λ) * ξ^-2 +
           -4B_W(κ, ϵ, Λ) * P_dξ_dκ(ξ, κ, ϵ, Λ) * ξ^-3 +
           6B_W(κ, ϵ, Λ) * P_dκ(ξ, κ, ϵ, Λ) * ξ^-4
end

function H(
    ξ,
    κ,
    ϵ,
    Λ::CGLParams;
    p = P(ξ, κ, ϵ, Λ),
    p_dϵ = P_dϵ(ξ, κ, ϵ, Λ),
    BW = B_W(κ, ϵ, Λ),
    BW_dϵ = B_W_dϵ(κ, ϵ, Λ),
)
    c_dϵ = _c_dϵ(κ, ϵ, Λ)

    return p * (-c_dϵ * BW + BW_dϵ * ξ^-2) + BW * p_dϵ * ξ^-2
end

function H_dξ(ξ, κ, ϵ, Λ::CGLParams)
    c_dϵ = _c_dϵ(κ, ϵ, Λ)

    return -c_dϵ * B_W(κ, ϵ, Λ) * P_dξ(ξ, κ, ϵ, Λ) +
           B_W_dϵ(κ, ϵ, Λ) * P_dξ(ξ, κ, ϵ, Λ) * ξ^-2 +
           -2B_W_dϵ(κ, ϵ, Λ) * P(ξ, κ, ϵ, Λ) * ξ^-3 +
           B_W(κ, ϵ, Λ) * P_dξ_dϵ(ξ, κ, ϵ, Λ) * ξ^-2 +
           -2B_W(κ, ϵ, Λ) * P_dϵ(ξ, κ, ϵ, Λ) * ξ^-3
end

function H_dξ_dξ(ξ, κ, ϵ, Λ::CGLParams)
    c_dϵ = _c_dϵ(κ, ϵ, Λ)

    return -c_dϵ * B_W(κ, ϵ, Λ) * P_dξ_dξ(ξ, κ, ϵ, Λ) +
           B_W_dϵ(κ, ϵ, Λ) * P_dξ_dξ(ξ, κ, ϵ, Λ) * ξ^-2 +
           -4B_W_dϵ(κ, ϵ, Λ) * P_dξ(ξ, κ, ϵ, Λ) * ξ^-3 +
           6B_W_dϵ(κ, ϵ, Λ) * P(ξ, κ, ϵ, Λ) * ξ^-4 +
           B_W(κ, ϵ, Λ) * P_dξ_dξ_dϵ(ξ, κ, ϵ, Λ) * ξ^-2 +
           -4B_W(κ, ϵ, Λ) * P_dξ_dϵ(ξ, κ, ϵ, Λ) * ξ^-3 +
           6B_W(κ, ϵ, Λ) * P_dϵ(ξ, κ, ϵ, Λ) * ξ^-4
end
