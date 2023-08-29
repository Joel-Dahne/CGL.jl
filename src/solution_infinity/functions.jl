export P, P_dξ, E, E_dξ, W

function P(ξ, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2

    return hypgeom_u(a, b, z)
end

# Used when automatic differentiation is wanted
function P_asym_approx(ξ, (p, κ)::Tuple{AbstractGLParams,Any})
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    b = oftype(a, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2

    return hypgeom_u_asym_approx(a, b, z)
end

function P_dξ(ξ, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    dzdξ = -im * κ / (1 - im * ϵ) * ξ

    return hypgeom_u_dz(a, b, z) * dzdξ
end

# Used when automatic differentiation is wanted
function P_dξ_asym_approx(ξ, (p, κ)::Tuple{AbstractGLParams,Any})
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    b = oftype(a, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    dzdξ = -im * κ / (1 - im * ϵ) * ξ

    return hypgeom_u_dz_asym_approx(a, b, z) * dzdξ
end

function P_dκ(ξ, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    a_dκ = -im * (ω / κ^2) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    z_dκ = -im / (1 - im * ϵ) * ξ^2 / 2

    return hypgeom_u_da(a, b, z) * a_dκ + hypgeom_u_dz(a, b, z) * z_dκ
end

function P_dξ_dκ(ξ, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    a_dκ = -im * (ω / κ^2) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    z_dξ = -im * κ / (1 - im * ϵ) * ξ
    z_dκ = -im / (1 - im * ϵ) * ξ^2 / 2
    z_dξ_dκ = -im / (1 - im * ϵ) * ξ

    return (hypgeom_u_dzda(a, b, z) * a_dκ + hypgeom_u_dz(a, b, z, 2) * z_dκ) * z_dξ +
           hypgeom_u_dz(a, b, z) * z_dξ_dκ
end

function E(ξ, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2

    return exp(z) * hypgeom_u(b - a, b, -z)
end

function E_dξ(ξ, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    z_dξ = -im * κ / (1 - im * ϵ) * ξ

    return exp(z) * (hypgeom_u(b - a, b, -z) - hypgeom_u_dz(b - a, b, -z)) * z_dξ
end

function E_dκ(ξ, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    a_dκ = -im * (ω / κ^2) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    z_dκ = -im / (1 - im * ϵ) * ξ^2 / 2

    return exp(z) * (
        z_dκ * hypgeom_u(b - a, b, -z) +
        (hypgeom_u_da(b - a, b, -z) * (-a_dκ) + hypgeom_u_dz(b - a, b, -z) * (z_dκ))
    )
end

function E_dξ_dκ(ξ, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    a_dκ = -im * (ω / κ^2) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    z_dξ = -im * κ / (1 - im * ϵ) * ξ
    z_dκ = -im / (1 - im * ϵ) * ξ^2 / 2
    z_dξ_dκ = -im / (1 - im * ϵ) * ξ

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

function W(ξ, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2

    # FIXME: The paper has ± in the first exp, what should we take?
    -im * κ / (1 - im * ϵ) * exp(im * (b - a) * π) * ξ * z^-b * exp(z)
end

function K(ξ, η, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    if η <= ξ
        return -1 / (1 - im * ϵ) * P(ξ, (p, κ)) * E(η, (p, κ)) / W(η, (p, κ))
    else
        return -1 / (1 - im * ϵ) * E(ξ, (p, κ)) * P(η, (p, κ)) / W(η, (p, κ))
    end
end
