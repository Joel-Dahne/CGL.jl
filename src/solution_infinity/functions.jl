export P, P_dξ, P_dκ, P_dξ_dκ, E, E_dξ, E_dκ, E_dξ_dκ, W, K, J_E, J_P, J_E_dκ, J_P_dκ

function _abc(κ, λ::AbstractGLParams{T}) where {T}
    (; d, ω, σ, ϵ) = λ

    a = (1 / σ + im * ω / κ) / 2
    b = convert(T, d) / 2
    c = -im * κ / 2(1 - im * ϵ)

    return a, b, c
end

function _abc(κ, λ::AbstractGLParams{Arb})
    (; d, ω, σ, ϵ) = λ

    a = (1 / σ + Acb(0, ω) / κ) / 2
    b = Acb(d // 2)
    c = Acb(0, -1) * κ / 2Acb(1, -ϵ)

    return a, b, c
end

function _abc_dκ(κ, λ::AbstractGLParams{T}) where {T}
    (; d, ω, σ, ϵ) = λ

    a = (1 / σ + im * ω / κ) / 2
    a_dκ = -im * (ω / κ^2) / 2
    b = convert(T, d) / 2
    c = -im * κ / (1 - im * ϵ) / 2
    c_dκ = -im / (1 - im * ϵ) / 2

    return a, a_dκ, b, c, c_dκ
end

function P(ξ, (λ, κ)::Tuple{AbstractGLParams{T},Any}) where {T}
    a, b, c = _abc(κ, λ)

    return hypgeom_u(a, b, c * ξ^2)
end

# Used when automatic differentiation is wanted
function P_asym_approx(ξ, (λ, κ)::Tuple{AbstractGLParams,Any})
    a, b, c = _abc(κ, λ)

    return hypgeom_u_asym_approx(a, b, c * ξ^2)
end

function P_dξ(ξ, (λ, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    a, b, c = _abc(κ, λ)

    return hypgeom_u_dz(a, b, c * ξ^2) * 2c * ξ
end

# Used when automatic differentiation is wanted
function P_dξ_asym_approx(ξ, (λ, κ)::Tuple{AbstractGLParams,Any})
    a, b, c = _abc(κ, λ)

    return hypgeom_u_dz_asym_approx(a, b, c * ξ^2) * 2c * ξ
end

function P_dκ(ξ, (λ, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    z = c * ξ^2
    z_dκ = c_dκ * ξ^2

    return hypgeom_u_da(a, b, z) * a_dκ + hypgeom_u_dz(a, b, z) * z_dκ
end

function P_dξ_dκ(ξ, (λ, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dκ = c_dκ * ξ^2
    z_dξ_dκ = 2c_dκ * ξ

    return (hypgeom_u_dzda(a, b, z) * a_dκ + hypgeom_u_dz(a, b, z, 2) * z_dκ) * z_dξ +
           hypgeom_u_dz(a, b, z) * z_dξ_dκ
end

function E(ξ, (λ, κ)::Tuple{AbstractGLParams{T},Any}) where {T}
    a, b, c = _abc(κ, λ)

    z = c * ξ^2

    return exp(z) * hypgeom_u(b - a, b, -z)
end

function E_dξ(ξ, (λ, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    a, b, c = _abc(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ

    return exp(z) * (hypgeom_u(b - a, b, -z) - hypgeom_u_dz(b - a, b, -z)) * z_dξ
end

function E_dκ(ξ, (λ, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    z = c * ξ^2
    z_dκ = c_dκ * ξ^2

    return exp(z) * (
        z_dκ * hypgeom_u(b - a, b, -z) +
        (hypgeom_u_da(b - a, b, -z) * (-a_dκ) + hypgeom_u_dz(b - a, b, -z) * (-z_dκ))
    )
end

function E_dξ_dκ(ξ, (λ, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ
    z_dκ = c_dκ * ξ^2
    z_dξ_dκ = 2c_dκ * ξ

    return exp(z) * (
        (hypgeom_u(b - a, b, -z) - hypgeom_u_dz(b - a, b, -z)) * z_dξ * z_dκ +
        (
            hypgeom_u_da(b - a, b, -z) * (-a_dκ) + hypgeom_u_dz(b - a, b, -z) * (-z_dκ) - (
                hypgeom_u_dzda(b - a, b, -z) * (-a_dκ) +
                hypgeom_u_dz(b - a, b, -z, 2) * (-z_dκ)
            )
        ) * z_dξ +
        (hypgeom_u(b - a, b, -z) - hypgeom_u_dz(b - a, b, -z)) * z_dξ_dκ
    )
end

function W(ξ, (λ, κ)::Tuple{AbstractGLParams{T},Any}) where {T}
    (; ϵ) = λ

    a, b, c = _abc(κ, λ)

    z = c * ξ^2

    # FIXME: The paper has ± in the first exp, what should we take?
    -im * κ / (1 - im * ϵ) * exp(im * (b - a) * π) * ξ * z^-b * exp(z)
end

function K(ξ, η, (λ, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    (; ϵ) = λ

    if η <= ξ
        return -1 / (1 - im * ϵ) * P(ξ, (λ, κ)) * E(η, (λ, κ)) / W(η, (λ, κ))
    else
        return -1 / (1 - im * ϵ) * E(ξ, (λ, κ)) * P(η, (λ, κ)) / W(η, (λ, κ))
    end
end

function J_P(ξ, (λ, κ)::Tuple{AbstractGLParams,Any})
    (; ϵ, δ) = λ

    return -(1 + im * δ) / (1 - im * ϵ) * P(ξ, (λ, κ)) * inv(W(ξ, (λ, κ)))
end

function J_E(ξ, (λ, κ)::Tuple{AbstractGLParams,Any})
    (; ϵ, δ) = λ

    return -(1 + im * δ) / (1 - im * ϵ) * E(ξ, (λ, κ)) * inv(W(ξ, (λ, κ)))
end

J_P_dκ(ξ, (λ, κ)::Tuple{AbstractGLParams,Any}) = J_P(ξ, (λ, ArbSeries((κ, 1))))[1]

J_E_dκ(ξ, (λ, κ)::Tuple{AbstractGLParams,Any}) = J_E(ξ, (λ, ArbSeries((κ, 1))))[1]
