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

struct FunctionEnclosures{T}
    P::T
    P_dξ::T
    P_dξ_dξ::T
    P_dκ::T
    P_dξ_dκ::T
    P_dϵ::T
    P_dξ_dϵ::T
    E::T
    E_dξ::T
    E_dκ::T
    E_dξ_dκ::T
    E_dϵ::T
    E_dξ_dϵ::T
    J_P::T
    J_P_dκ::T
    J_P_dϵ::T
    J_E::T
    J_E_dκ::T
    J_E_dϵ::T
    D::T
    D_dξ::T
    H::T
    H_dξ::T

    FunctionEnclosures(
        ξ::Arb,
        κ::Arb,
        ϵ::Arb,
        λ::CGLParams{Arb};
        include_dκ::Bool = false,
        include_dϵ::Bool = false,
    ) = FunctionEnclosures{Acb}(ξ, κ, ϵ, λ; include_dκ, include_dϵ)

    FunctionEnclosures(
        ξ::S,
        κ::S,
        ϵ::S,
        λ::CGLParams{S};
        include_dκ::Bool = false,
        include_dϵ::Bool = false,
    ) where {S} = FunctionEnclosures{complex(S)}(ξ, κ, ϵ, λ; include_dκ, include_dϵ)

    function FunctionEnclosures{T}(
        ξ::S,
        κ::S,
        ϵ::S,
        λ::CGLParams{S};
        include_dκ::Bool = false,
        include_dϵ::Bool = false,
    ) where {S,T}
        return new(
            P(ξ, κ, ϵ, λ),
            P_dξ(ξ, κ, ϵ, λ),
            P_dξ_dξ(ξ, κ, ϵ, λ),
            include_dκ ? P_dκ(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dκ ? P_dξ_dκ(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dϵ ? P_dϵ(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dϵ ? P_dξ_dϵ(ξ, κ, ϵ, λ) : indeterminate(T),
            E(ξ, κ, ϵ, λ),
            E_dξ(ξ, κ, ϵ, λ),
            include_dκ ? E_dκ(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dκ ? E_dξ_dκ(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dϵ ? E_dϵ(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dϵ ? E_dξ_dϵ(ξ, κ, ϵ, λ) : indeterminate(T),
            J_P(ξ, κ, ϵ, λ),
            include_dκ ? J_P_dκ(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dϵ ? J_P_dϵ(ξ, κ, ϵ, λ) : indeterminate(T),
            J_E(ξ, κ, ϵ, λ),
            include_dκ ? J_E_dκ(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dϵ ? J_E_dϵ(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dκ ? D(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dκ ? D_dξ(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dϵ ? H(ξ, κ, ϵ, λ) : indeterminate(T),
            include_dϵ ? H_dξ(ξ, κ, ϵ, λ) : indeterminate(T),
        )
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

J_P(ξ, κ, ϵ, λ::CGLParams) = (1 + im * λ.δ) / (1 - im * ϵ) * P(ξ, κ, ϵ, λ) / W(ξ, κ, ϵ, λ)

J_E(ξ, κ, ϵ, λ::CGLParams) = (1 + im * λ.δ) / (1 - im * ϵ) * E(ξ, κ, ϵ, λ) / W(ξ, κ, ϵ, λ)

# IMPROVE: These always return Acb even if the input is for example
# Float64.

J_P_dξ(ξ, κ, ϵ, λ::CGLParams) = J_P(ArbSeries((ξ, 1)), κ, ϵ, λ)[1]

J_E_dξ(ξ, κ, ϵ, λ::CGLParams) = J_E(ArbSeries((ξ, 1)), κ, ϵ, λ)[1]

J_P_dξ_dξ(ξ, κ, ϵ, λ::CGLParams) = 2J_P(ArbSeries((ξ, 1), degree = 2), κ, ϵ, λ)[2]

J_E_dξ_dξ(ξ, κ, ϵ, λ::CGLParams) = 2J_E(ArbSeries((ξ, 1), degree = 2), κ, ϵ, λ)[2]

J_P_dκ(ξ, κ, ϵ, λ::CGLParams) = J_P(ξ, ArbSeries((κ, 1)), ϵ, λ)[1]

J_E_dκ(ξ, κ, ϵ, λ::CGLParams) = J_E(ξ, ArbSeries((κ, 1)), ϵ, λ)[1]

J_P_dϵ(ξ, κ, ϵ, λ::CGLParams) = J_P(ξ, κ, ArbSeries((ϵ, 1)), λ)[1]

J_E_dϵ(ξ, κ, ϵ, λ::CGLParams) = J_E(ξ, κ, ArbSeries((ϵ, 1)), λ)[1]

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

function D(ξ, κ, ϵ, λ::CGLParams)
    _, _, _, _, c_dκ = _abc_dκ(κ, ϵ, λ)

    return -c_dκ * B_W(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) +
           B_W_dκ(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) * ξ^-2 +
           B_W(κ, ϵ, λ) * P_dκ(ξ, κ, ϵ, λ) * ξ^-2
end

function D_dξ(ξ, κ, ϵ, λ::CGLParams)
    _, _, _, _, c_dκ = _abc_dκ(κ, ϵ, λ)

    return -c_dκ * B_W(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) +
           B_W_dκ(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) * ξ^-2 +
           -2B_W_dκ(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) * ξ^-3 +
           B_W(κ, ϵ, λ) * P_dξ_dκ(ξ, κ, ϵ, λ) * ξ^-2 +
           -2B_W(κ, ϵ, λ) * P_dκ(ξ, κ, ϵ, λ) * ξ^-3
end

function D_dξ_dξ(ξ, κ, ϵ, λ::CGLParams)
    _, _, _, _, c_dκ = _abc_dκ(κ, ϵ, λ)

    return -c_dκ * B_W(κ, ϵ, λ) * P_dξ_dξ(ξ, κ, ϵ, λ) +
           B_W_dκ(κ, ϵ, λ) * P_dξ_dξ(ξ, κ, ϵ, λ) * ξ^-2 +
           -4B_W_dκ(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) * ξ^-3 +
           6B_W_dκ(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) * ξ^-4 +
           B_W(κ, ϵ, λ) * P_dξ_dξ_dκ(ξ, κ, ϵ, λ) * ξ^-2 +
           -4B_W(κ, ϵ, λ) * P_dξ_dκ(ξ, κ, ϵ, λ) * ξ^-3 +
           6B_W(κ, ϵ, λ) * P_dκ(ξ, κ, ϵ, λ) * ξ^-4
end

#
function H(ξ, κ, ϵ, λ::CGLParams)
    _, _, _, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    return -c_dϵ * B_W(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) +
           B_W_dϵ(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) * ξ^-2 +
           B_W(κ, ϵ, λ) * P_dϵ(ξ, κ, ϵ, λ) * ξ^-2
end

function H_dξ(ξ, κ, ϵ, λ::CGLParams)
    _, _, _, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    return -c_dϵ * B_W(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) +
           B_W_dϵ(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) * ξ^-2 +
           -2B_W_dϵ(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) * ξ^-3 +
           B_W(κ, ϵ, λ) * P_dξ_dϵ(ξ, κ, ϵ, λ) * ξ^-2 +
           -2B_W(κ, ϵ, λ) * P_dϵ(ξ, κ, ϵ, λ) * ξ^-3
end

function H_dξ_dξ(ξ, κ, ϵ, λ::CGLParams)
    _, _, _, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    return -c_dϵ * B_W(κ, ϵ, λ) * P_dξ_dξ(ξ, κ, ϵ, λ) +
           B_W_dϵ(κ, ϵ, λ) * P_dξ_dξ(ξ, κ, ϵ, λ) * ξ^-2 +
           -4B_W_dϵ(κ, ϵ, λ) * P_dξ(ξ, κ, ϵ, λ) * ξ^-3 +
           6B_W_dϵ(κ, ϵ, λ) * P(ξ, κ, ϵ, λ) * ξ^-4 +
           B_W(κ, ϵ, λ) * P_dξ_dξ_dϵ(ξ, κ, ϵ, λ) * ξ^-2 +
           -4B_W(κ, ϵ, λ) * P_dξ_dϵ(ξ, κ, ϵ, λ) * ξ^-3 +
           6B_W(κ, ϵ, λ) * P_dϵ(ξ, κ, ϵ, λ) * ξ^-4
end
