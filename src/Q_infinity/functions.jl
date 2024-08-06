export P, P_dξ, P_dξ_dξ, P_dξ_dξ_dξ
export P_dκ, P_dξ_dκ, P_dξ_dξ_dκ
export P_dϵ, P_dξ_dϵ, P_dξ_dξ_dϵ

export E, E_dξ, E_dξ_dξ, E_dξ_dξ_dξ
export E_dκ, E_dξ_dκ
export E_dϵ, E_dξ_dϵ

export W

export J_E, J_E_dξ, J_E_dξ_dξ, J_E_dκ, J_E_dϵ
export J_P, J_P_dξ, J_P_dξ_dξ, J_P_dκ, J_P_dϵ
export D, D_dξ, D_dξ_dξ, H, H_dξ, H_dξ_dξ

struct FunctionEnclosures
    P::Acb
    P_dξ::Acb
    P_dξ_dξ::Acb
    P_dκ::Acb
    P_dξ_dκ::Acb
    P_dϵ::Acb
    P_dξ_dϵ::Acb
    E::Acb
    E_dξ::Acb
    E_dκ::Acb
    E_dξ_dκ::Acb
    E_dϵ::Acb
    E_dξ_dϵ::Acb
    J_P::Acb
    J_P_dκ::Acb
    J_P_dϵ::Acb
    J_E::Acb
    J_E_dκ::Acb
    J_E_dϵ::Acb
    D::Acb
    D_dξ::Acb
    H::Acb
    H_dξ::Acb

    function FunctionEnclosures(
        ξ::Arb,
        κ::Arb,
        ϵ::Arb,
        λ::CGLParams{Arb};
        include_dκ::Bool = false,
        include_dϵ::Bool = false,
    )
        ξ_series = ArbSeries((ξ, 1))
        BW = B_W(κ, ϵ, λ)

        P_series = P(ArbSeries(ξ_series, degree = 2), κ, ϵ, λ)
        p = P_series[0]
        p_dξ = P_series[1]
        p_dξ_dξ = 2P_series[2]

        E_series = E(ξ_series, κ, ϵ, λ)
        e = E_series[0]
        e_dξ = E_series[1]

        if include_dκ
            BW_dκ = B_W_dκ(κ, ϵ, λ)

            P_dκ_series = P_dκ(ξ_series, κ, ϵ, λ)
            p_dκ = P_dκ_series[0]
            p_dξ_dκ = P_dκ_series[1]

            E_dκ_series = E_dκ(ξ_series, κ, ϵ, λ)
            e_dκ = E_dκ_series[0]
            e_dξ_dκ = E_dκ_series[1]

            j_e_dκ = J_E_dκ(ξ, κ, ϵ, λ; e, e_dκ, BW, BW_dκ)

            D_series = D(ξ_series, κ, ϵ, λ; p = P_series, p_dκ = P_dκ_series, BW, BW_dκ)
            d = D_series[0]
            d_dξ = D_series[1]
        else
            p_dκ = indeterminate(Acb)
            p_dξ_dκ = indeterminate(Acb)

            e_dκ = indeterminate(Acb)
            e_dξ_dκ = indeterminate(Acb)

            j_e_dκ = indeterminate(Acb)

            d = indeterminate(Acb)
            d_dξ = indeterminate(Acb)
        end

        if include_dϵ
            BW_dϵ = B_W_dϵ(κ, ϵ, λ)

            P_dϵ_series = P_dϵ(ξ_series, κ, ϵ, λ)
            p_dϵ = P_dϵ_series[0]
            p_dξ_dϵ = P_dϵ_series[1]

            E_dϵ_series = E_dϵ(ξ_series, κ, ϵ, λ)
            e_dϵ = E_dϵ_series[0]
            e_dξ_dϵ = E_dϵ_series[1]

            j_e_dϵ = J_E_dϵ(ξ, κ, ϵ, λ; e, e_dϵ, BW, BW_dϵ)

            H_series = H(ξ_series, κ, ϵ, λ; p = P_series, p_dϵ = P_dϵ_series, BW, BW_dϵ)
            h = H_series[0]
            h_dξ = H_series[1]
        else
            p_dϵ = indeterminate(Acb)
            p_dξ_dϵ = indeterminate(Acb)

            e_dϵ = indeterminate(Acb)
            e_dξ_dϵ = indeterminate(Acb)

            j_e_dϵ = indeterminate(Acb)

            h = indeterminate(Acb)
            h_dξ = indeterminate(Acb)
        end

        j_p_dκ =
            F = new(
                p,
                p_dξ,
                p_dξ_dξ,
                p_dκ,
                p_dξ_dκ,
                p_dϵ,
                p_dξ_dϵ,
                e,
                e_dξ,
                e_dκ,
                e_dξ_dκ,
                e_dϵ,
                e_dξ_dϵ,
                J_P(ξ, κ, ϵ, λ; p),
                include_dκ ? J_P_dκ(ξ, κ, ϵ, λ; d) : indeterminate(Acb),
                include_dϵ ? J_P_dϵ(ξ, κ, ϵ, λ; h) : indeterminate(Acb),
                J_E(ξ, κ, ϵ, λ; e),
                j_e_dκ,
                j_e_dϵ,
                d,
                d_dξ,
                h,
                h_dξ,
            )

        return F
    end
end

function P(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = c * ξ^2

    return U(a, b, z)
end

function P_dξ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ

    return U_dz(a, b, z) * z_dξ
end

function P_dξ_dξ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c

    return U_dz(a, b, z, 2) * z_dξ^2 + U_dz(a, b, z) * z_dξ_dξ
end

function P_dξ_dξ_dξ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c
    # z_dξ_dξ_dξ = 0

    return U_dz(a, b, z, 3) * z_dξ^3 + U_dz(a, b, z, 2) * 3z_dξ * z_dξ_dξ
end

function P_dκ(ξ, κ, ϵ, λ::CGLParams)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    z = c * ξ^2
    z_dκ = c_dκ * ξ^2

    return U_da(a, b, z) * a_dκ + U_dz(a, b, z) * z_dκ
end

function P_dξ_dκ(ξ, κ, ϵ, λ::CGLParams)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dκ = c_dκ * ξ^2
    z_dξ_dκ = 2c_dκ * ξ

    return (U_dzda(a, b, z) * a_dκ + U_dz(a, b, z, 2) * z_dκ) * z_dξ +
           U_dz(a, b, z) * z_dξ_dκ
end

function P_dξ_dξ_dκ(ξ, κ, ϵ, λ::CGLParams)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

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

function P_dϵ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    z = c * ξ^2
    z_dϵ = c_dϵ * ξ^2

    return U_dz(a, b, z) * z_dϵ
end

function P_dξ_dϵ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dϵ = c_dϵ * ξ^2
    z_dξ_dϵ = 2c_dϵ * ξ

    return U_dz(a, b, z, 2) * z_dϵ * z_dξ + U_dz(a, b, z) * z_dξ_dϵ
end

function P_dξ_dξ_dϵ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

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

function E(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = c * ξ^2

    return exp(z) * U(b - a, b, -z)
end

function E_dξ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ

    return exp(z) * (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ
end

function E_dξ_dξ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c

    return exp(z) * (
        U(b - a, b, -z) * (z_dξ^2 + z_dξ_dξ) - U_dz(b - a, b, -z) * (2z_dξ^2 + z_dξ_dξ) +
        U_dz(b - a, b, -z, 2) * z_dξ^2
    )
end

function E_dξ_dξ_dξ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

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

function E_dκ(ξ, κ, ϵ, λ::CGLParams)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

    z = c * ξ^2
    z_dκ = c_dκ * ξ^2

    return exp(z) * (
        z_dκ * U(b - a, b, -z) +
        (U_da(b - a, b, -z) * (-a_dκ) + U_dz(b - a, b, -z) * (-z_dκ))
    )
end

function E_dξ_dκ(ξ, κ, ϵ, λ::CGLParams)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, ϵ, λ)

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

function E_dϵ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    z = c * ξ^2
    z_dϵ = c_dϵ * ξ^2

    return exp(z) * (z_dϵ * U(b - a, b, -z) + U_dz(b - a, b, -z) * (-z_dϵ))
end

function E_dξ_dϵ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

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

function W(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = c * ξ^2

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 2c * exp(sgn * im * (b - a) * π) * ξ * z^-b * exp(z)
end

function J_P(ξ, κ, ϵ, λ::CGLParams; p = P(ξ, κ, ϵ, λ))
    c = _c(κ, ϵ, λ)

    return B_W(κ, ϵ, λ) * p * exp(-c * ξ^2) * ξ^(λ.d - 1)
end

function J_E(ξ, κ, ϵ, λ::CGLParams; e = E(ξ, κ, ϵ, λ))
    c = _c(κ, ϵ, λ)

    return B_W(κ, ϵ, λ) * e * exp(-c * ξ^2) * ξ^(λ.d - 1)
end

# These four are only used for testing and are not performance critical
J_P_dξ(ξ, κ, ϵ, λ::CGLParams) = J_P(ArbSeries((ξ, 1)), κ, ϵ, λ)[1]

J_E_dξ(ξ, κ, ϵ, λ::CGLParams) = J_E(ArbSeries((ξ, 1)), κ, ϵ, λ)[1]

J_P_dξ_dξ(ξ, κ, ϵ, λ::CGLParams) = 2J_P(ArbSeries((ξ, 1), degree = 2), κ, ϵ, λ)[2]

J_E_dξ_dξ(ξ, κ, ϵ, λ::CGLParams) = 2J_E(ArbSeries((ξ, 1), degree = 2), κ, ϵ, λ)[2]

function J_P_dκ(ξ, κ, ϵ, λ::CGLParams; d = D(ξ, κ, ϵ, λ))
    c = _c(κ, ϵ, λ)

    return d * exp(-c * ξ^2) * ξ^(λ.d + 1)
end

function J_E_dκ(
    ξ,
    κ,
    ϵ,
    λ::CGLParams;
    e = E(ξ, κ, ϵ, λ),
    e_dκ = E_dκ(ξ, κ, ϵ, λ),
    BW = B_W(κ, ϵ, λ),
    BW_dκ = B_W_dκ(κ, ϵ, λ),
)
    c, c_dκ = _c(κ, ϵ, λ), _c_dκ(κ, ϵ, λ)

    return (BW_dκ * e * ξ^-2 + BW * e_dκ * ξ^-2 - BW * e * c_dκ) *
           exp(-c * ξ^2) *
           ξ^(λ.d + 1)
end

function J_P_dϵ(ξ, κ, ϵ, λ::CGLParams; h = H(ξ, κ, ϵ, λ))
    c = _c(κ, ϵ, λ)

    return h * exp(-c * ξ^2) * ξ^(λ.d + 1)
end

function J_E_dϵ(
    ξ,
    κ,
    ϵ,
    λ::CGLParams;
    e = E(ξ, κ, ϵ, λ),
    e_dϵ = E_dϵ(ξ, κ, ϵ, λ),
    BW = B_W(κ, ϵ, λ),
    BW_dϵ = B_W_dϵ(κ, ϵ, λ),
)
    c, c_dϵ = _c(κ, ϵ, λ), _c_dϵ(κ, ϵ, λ)

    return (BW_dϵ * e * ξ^-2 + BW * e_dϵ * ξ^-2 - BW * e * c_dϵ) *
           exp(-c * ξ^2) *
           ξ^(λ.d + 1)
end

J_P_dξ(ξ::Float64, κ::Float64, ϵ::Float64, λ::CGLParams{Float64}) =
    ComplexF64(J_P(ArbSeries((ξ, 1)), κ, ϵ, λ)[1])
J_E_dξ(ξ::Float64, κ::Float64, ϵ::Float64, λ::CGLParams{Float64}) =
    ComplexF64(J_E(ArbSeries((ξ, 1)), κ, ϵ, λ)[1])
J_P_dξ_dξ(ξ::Float64, κ::Float64, ϵ::Float64, λ::CGLParams{Float64}) =
    ComplexF64(2J_P(ArbSeries((ξ, 1), degree = 2), κ, ϵ, λ)[2])
J_E_dξ_dξ(ξ::Float64, κ::Float64, ϵ::Float64, λ::CGLParams{Float64}) =
    ComplexF64(2J_E(ArbSeries((ξ, 1), degree = 2), κ, ϵ, λ)[2])
J_P_dκ(ξ::Float64, κ::Float64, ϵ::Float64, λ::CGLParams{Float64}) =
    ComplexF64(J_P(ξ, ArbSeries((κ, 1)), ϵ, λ)[1])
J_E_dκ(ξ::Float64, κ::Float64, ϵ::Float64, λ::CGLParams{Float64}) =
    ComplexF64(J_E(ξ, ArbSeries((κ, 1)), ϵ, λ)[1])
J_P_dϵ(ξ::Float64, κ::Float64, ϵ::Float64, λ::CGLParams{Float64}) =
    ComplexF64(J_P(ξ, κ, ArbSeries((ϵ, 1)), λ)[1])
J_E_dϵ(ξ::Float64, κ::Float64, ϵ::Float64, λ::CGLParams{Float64}) =
    ComplexF64(J_E(ξ, κ, ArbSeries((ϵ, 1)), λ)[1])

function D(
    ξ,
    κ,
    ϵ,
    λ::CGLParams;
    p = P(ξ, κ, ϵ, λ),
    p_dκ = P_dκ(ξ, κ, ϵ, λ),
    BW = B_W(κ, ϵ, λ),
    BW_dκ = B_W_dκ(κ, ϵ, λ),
)
    c_dκ = _c_dκ(κ, ϵ, λ)

    return p * (-c_dκ * BW + BW_dκ * ξ^-2) + BW * p_dκ * ξ^-2
end

function D_dξ(ξ, κ, ϵ, λ::CGLParams)
    c_dκ = _c_dκ(κ, ϵ, λ)

    return -c_dκ * B_W(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) +
           B_W_dκ(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) * ξ^-2 +
           -2B_W_dκ(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) * ξ^-3 +
           B_W(κ, ϵ, λ) * P_dξ_dκ(ξ, κ, ϵ, λ) * ξ^-2 +
           -2B_W(κ, ϵ, λ) * P_dκ(ξ, κ, ϵ, λ) * ξ^-3
end

function D_dξ_dξ(ξ, κ, ϵ, λ::CGLParams)
    c_dκ = _c_dκ(κ, ϵ, λ)

    return -c_dκ * B_W(κ, ϵ, λ) * P_dξ_dξ(ξ, κ, ϵ, λ) +
           B_W_dκ(κ, ϵ, λ) * P_dξ_dξ(ξ, κ, ϵ, λ) * ξ^-2 +
           -4B_W_dκ(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) * ξ^-3 +
           6B_W_dκ(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) * ξ^-4 +
           B_W(κ, ϵ, λ) * P_dξ_dξ_dκ(ξ, κ, ϵ, λ) * ξ^-2 +
           -4B_W(κ, ϵ, λ) * P_dξ_dκ(ξ, κ, ϵ, λ) * ξ^-3 +
           6B_W(κ, ϵ, λ) * P_dκ(ξ, κ, ϵ, λ) * ξ^-4
end

function H(
    ξ,
    κ,
    ϵ,
    λ::CGLParams;
    p = P(ξ, κ, ϵ, λ),
    p_dϵ = P_dϵ(ξ, κ, ϵ, λ),
    BW = B_W(κ, ϵ, λ),
    BW_dϵ = B_W_dϵ(κ, ϵ, λ),
)
    c_dϵ = _c_dϵ(κ, ϵ, λ)

    return p * (-c_dϵ * BW + BW_dϵ * ξ^-2) + BW * p_dϵ * ξ^-2
end

function H_dξ(ξ, κ, ϵ, λ::CGLParams)
    c_dϵ = _c_dϵ(κ, ϵ, λ)

    return -c_dϵ * B_W(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) +
           B_W_dϵ(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) * ξ^-2 +
           -2B_W_dϵ(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) * ξ^-3 +
           B_W(κ, ϵ, λ) * P_dξ_dϵ(ξ, κ, ϵ, λ) * ξ^-2 +
           -2B_W(κ, ϵ, λ) * P_dϵ(ξ, κ, ϵ, λ) * ξ^-3
end

function H_dξ_dξ(ξ, κ, ϵ, λ::CGLParams)
    c_dϵ = _c_dϵ(κ, ϵ, λ)

    return -c_dϵ * B_W(κ, ϵ, λ) * P_dξ_dξ(ξ, κ, ϵ, λ) +
           B_W_dϵ(κ, ϵ, λ) * P_dξ_dξ(ξ, κ, ϵ, λ) * ξ^-2 +
           -4B_W_dϵ(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) * ξ^-3 +
           6B_W_dϵ(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) * ξ^-4 +
           B_W(κ, ϵ, λ) * P_dξ_dξ_dϵ(ξ, κ, ϵ, λ) * ξ^-2 +
           -4B_W(κ, ϵ, λ) * P_dξ_dϵ(ξ, κ, ϵ, λ) * ξ^-3 +
           6B_W(κ, ϵ, λ) * P_dϵ(ξ, κ, ϵ, λ) * ξ^-4
end
