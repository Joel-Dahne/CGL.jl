export P,
    P_dξ,
    P_dξ_dξ,
    P_dξ_dξ_dξ,
    P_dκ,
    P_dξ_dκ,
    P_dξ_dξ_dκ,
    E,
    E_dξ,
    E_dξ_dξ,
    E_dξ_dξ_dξ,
    E_dκ,
    E_dξ_dκ,
    W,
    J_E,
    J_P,
    J_E_dξ,
    J_P_dξ,
    J_E_dξ_dξ,
    J_P_dξ_dξ,
    J_E_dκ,
    J_P_dκ,
    D,
    D_dξ,
    D_dξ_dξ

function _abc(κ, λ::CGLParams{T}) where {T}
    (; d, ω, σ, ϵ) = λ

    a = (1 / σ + im * ω / κ) / 2
    b = convert(T, d) / 2
    c = -im * κ / 2(1 - im * ϵ)

    return a, b, c
end

function _abc(κ::Arb, λ::CGLParams{Arb})
    (; d, ω, σ, ϵ) = λ

    # Acb(1 / σ, ω / κ) / 2
    a = let a = Acb()
        Arblib.inv!(Arblib.realref(a), σ)
        Arblib.div!(Arblib.imagref(a), ω, κ)
        Arblib.mul_2exp!(a, a, -1)
    end
    b = Acb(d // 2)
    # c = κ / 2Acb(ϵ, 1)
    c = let c = Acb(κ)
        Arblib.mul_2exp!(c, c, -1)
        Arblib.div!(c, c, Acb(ϵ, 1))
    end

    return a, b, c
end

function _abc(κ::ArbSeries, λ::CGLParams{Arb})
    (; d, ω, σ, ϵ) = λ

    a = (1 / σ + im * ω / κ) / 2
    b = Acb(d // 2)
    c = κ / 2Acb(ϵ, 1)

    return a, b, c
end

function _abc_dκ(κ, λ::CGLParams{T}) where {T}
    (; d, ω, σ, ϵ) = λ

    a, b, c = _abc(κ, λ)

    if T == Arb
        a_dκ = Acb(0, -1) * (ω / κ^2) / 2
        c_dκ = Acb(0, -1) / Acb(1, -ϵ) / 2
    else
        a_dκ = -im * (ω / κ^2) / 2
        c_dκ = -im / (1 - im * ϵ) / 2
    end

    return a, a_dκ, b, c, c_dκ
end

function P(ξ, (λ, κ)::Tuple{CGLParams{T},Any}) where {T}
    a, b, c = _abc(κ, λ)

    z = c * ξ^2

    return U(a, b, z)
end

# Used when automatic differentiation is wanted
function P_asym_approx(ξ, (λ, κ)::Tuple{CGLParams,Any})
    a, b, c = _abc(κ, λ)

    z = c * ξ^2

    return U_asym_approx(a, b, z)
end

function P_dξ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, b, c = _abc(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ

    return U_dz(a, b, z) * z_dξ
end

# Used when automatic differentiation is wanted
function P_dξ_asym_approx(ξ, (λ, κ)::Tuple{CGLParams,Any})
    a, b, c = _abc(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ

    return U_dz_asym_approx(a, b, z) * z_dξ
end

function P_dξ_dξ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, b, c = _abc(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c

    return U_dz(a, b, z, 2) * z_dξ^2 + U_dz(a, b, z) * z_dξ_dξ
end

function P_dξ_dξ_dξ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, b, c = _abc(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c
    # z_dξ_dξ_dξ = 0

    return U_dz(a, b, z, 3) * z_dξ^3 + U_dz(a, b, z, 2) * 3z_dξ * z_dξ_dξ
end

function P_dκ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    z = c * ξ^2
    z_dκ = c_dκ * ξ^2

    return U_da(a, b, z) * a_dκ + U_dz(a, b, z) * z_dκ
end

function P_dξ_dκ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dκ = c_dκ * ξ^2
    z_dξ_dκ = 2c_dκ * ξ

    return (U_dzda(a, b, z) * a_dκ + U_dz(a, b, z, 2) * z_dκ) * z_dξ +
           U_dz(a, b, z) * z_dξ_dκ
end

function P_dξ_dξ_dκ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

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

function E(ξ, (λ, κ)::Tuple{CGLParams{T},Any}) where {T}
    a, b, c = _abc(κ, λ)

    z = c * ξ^2

    return exp(z) * U(b - a, b, -z)
end

function E_dξ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, b, c = _abc(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ

    return exp(z) * (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ
end

function E_dξ_dξ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, b, c = _abc(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dξ_dξ = 2c

    return exp(z) * (
        U(b - a, b, -z) * (z_dξ^2 + z_dξ_dξ) - U_dz(b - a, b, -z) * (2z_dξ^2 + z_dξ_dξ) +
        U_dz(b - a, b, -z, 2) * z_dξ^2
    )
end

function E_dξ_dξ_dξ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, b, c = _abc(κ, λ)

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

function E_dκ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    z = c * ξ^2
    z_dκ = c_dκ * ξ^2

    return exp(z) * (
        z_dκ * U(b - a, b, -z) +
        (U_da(b - a, b, -z) * (-a_dκ) + U_dz(b - a, b, -z) * (-z_dκ))
    )
end

function E_dξ_dκ(ξ, (λ, κ)::Tuple{CGLParams{T},T}) where {T}
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

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

function W(ξ, (λ, κ)::Tuple{CGLParams{T},Any}) where {T}
    (; ϵ) = λ

    a, b, c = _abc(κ, λ)

    z = c * ξ^2

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 2c * exp(sgn * im * (b - a) * π) * ξ * z^-b * exp(z)
end

function B_W(κ, λ::CGLParams{T}) where {T}
    (; δ) = λ

    a, b, c = _abc(κ, λ)

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    if T == Arb
        return -Acb(1, δ) / Acb(0, κ) * exp(-sgn * im * (b - a) * π) * c^b
    else
        return -(1 + im * δ) / (im * κ) * exp(-sgn * im * (b - a) * π) * c^b
    end
end

function B_W_dκ(κ::Arb, λ::CGLParams)
    (; δ) = λ

    κ_series = ArbSeries((κ, 1))

    a, b, c = _abc(κ_series, λ)

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    res = -Acb(1, δ) / (im * κ_series) * exp(-sgn * im * (b - a) * π) * c^b

    return res[1]
end

B_W_dκ(κ, λ) = ForwardDiff.derivative(κ -> B_W(κ, λ), κ)

function J_P(ξ, (λ, κ)::Tuple{CGLParams,Any})
    (; ϵ, δ) = λ

    return (1 + im * δ) / (1 - im * ϵ) * P(ξ, (λ, κ)) / W(ξ, (λ, κ))
end

function J_E(ξ, (λ, κ)::Tuple{CGLParams,Any})
    (; ϵ, δ) = λ

    return (1 + im * δ) / (1 - im * ϵ) * E(ξ, (λ, κ)) / W(ξ, (λ, κ))
end

# IMPROVE: These always return Acb even if the input is for example
# Float64.

J_P_dξ(ξ, (λ, κ)::Tuple{CGLParams,Any}) = J_P(ArbSeries((ξ, 1)), (λ, κ))[1]

J_E_dξ(ξ, (λ, κ)::Tuple{CGLParams,Any}) = J_E(ArbSeries((ξ, 1)), (λ, κ))[1]

J_P_dξ_dξ(ξ, (λ, κ)::Tuple{CGLParams,Any}) = 2J_P(ArbSeries((ξ, 1), degree = 2), (λ, κ))[2]

J_E_dξ_dξ(ξ, (λ, κ)::Tuple{CGLParams,Any}) = 2J_E(ArbSeries((ξ, 1), degree = 2), (λ, κ))[2]

J_P_dκ(ξ, (λ, κ)::Tuple{CGLParams,Any}) = J_P(ξ, (λ, ArbSeries((κ, 1))))[1]

J_E_dκ(ξ, (λ, κ)::Tuple{CGLParams,Any}) = J_E(ξ, (λ, ArbSeries((κ, 1))))[1]

J_P_dξ(ξ::Float64, (λ, κ)::Tuple{CGLParams{Float64},Float64}) =
    ComplexF64(J_P(ArbSeries((ξ, 1)), (λ, κ))[1])
J_E_dξ(ξ::Float64, (λ, κ)::Tuple{CGLParams{Float64},Float64}) =
    ComplexF64(J_E(ArbSeries((ξ, 1)), (λ, κ))[1])
J_P_dξ_dξ(ξ::Float64, (λ, κ)::Tuple{CGLParams{Float64},Float64}) =
    ComplexF64(2J_P(ArbSeries((ξ, 1), degree = 2), (λ, κ))[2])
J_E_dξ_dξ(ξ::Float64, (λ, κ)::Tuple{CGLParams{Float64},Float64}) =
    ComplexF64(2J_E(ArbSeries((ξ, 1), degree = 2), (λ, κ))[2])
J_P_dκ(ξ::Float64, (λ, κ)::Tuple{CGLParams{Float64},Float64}) =
    ComplexF64(J_P(ξ, (λ, ArbSeries((κ, 1))))[1])
J_E_dκ(ξ::Float64, (λ, κ)::Tuple{CGLParams{Float64},Float64}) =
    ComplexF64(J_E(ξ, (λ, ArbSeries((κ, 1))))[1])

function D(ξ, (λ, κ)::Tuple{CGLParams,Any})
    _, _, _, _, c_dκ = _abc_dκ(κ, λ)

    return -c_dκ * B_W(κ, λ) * P(ξ, (λ, κ)) +
           B_W_dκ(κ, λ) * P(ξ, (λ, κ)) * ξ^-2 +
           B_W(κ, λ) * P_dκ(ξ, (λ, κ)) * ξ^-2
end

function D_dξ(ξ, (λ, κ)::Tuple{CGLParams,Any})
    _, _, _, _, c_dκ = _abc_dκ(κ, λ)

    return -c_dκ * B_W(κ, λ) * P_dξ(ξ, (λ, κ)) +
           B_W_dκ(κ, λ) * P_dξ(ξ, (λ, κ)) * ξ^-2 +
           -2B_W_dκ(κ, λ) * P(ξ, (λ, κ)) * ξ^-3 +
           B_W(κ, λ) * P_dξ_dκ(ξ, (λ, κ)) * ξ^-2 +
           -2B_W(κ, λ) * P_dκ(ξ, (λ, κ)) * ξ^-3
end

function D_dξ_dξ(ξ, (λ, κ)::Tuple{CGLParams,Any})
    _, _, _, _, c_dκ = _abc_dκ(κ, λ)

    return -c_dκ * B_W(κ, λ) * P_dξ_dξ(ξ, (λ, κ)) +
           B_W_dκ(κ, λ) * P_dξ_dξ(ξ, (λ, κ)) * ξ^-2 +
           -4B_W_dκ(κ, λ) * P_dξ(ξ, (λ, κ)) * ξ^-3 +
           6B_W_dκ(κ, λ) * P(ξ, (λ, κ)) * ξ^-4 +
           B_W(κ, λ) * P_dξ_dξ_dκ(ξ, (λ, κ)) * ξ^-2 +
           -4B_W(κ, λ) * P_dξ_dκ(ξ, (λ, κ)) * ξ^-3 +
           6B_W(κ, λ) * P_dκ(ξ, (λ, κ)) * ξ^-4
end
