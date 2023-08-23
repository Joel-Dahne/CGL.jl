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
    dadκ = -im * (ω / κ^2) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    dzdκ = -im / (1 - im * ϵ) * ξ^2 / 2

    return hypgeom_u_da(a, b, z) * dadκ + hypgeom_u_dz(a, b, z) * dzdκ
end

function P_dξdκ(ξ, (p, κ)::Tuple{AbstractGLParams{T},T}) where {T}
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = (1 / σ + im * ω / κ) / 2
    dadκ = -im * (ω / κ^2) / 2
    b = convert(T, d) / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    dzdξ = -im * κ / (1 - im * ϵ) * ξ
    dzdκ = -im / (1 - im * ϵ) * ξ^2 / 2
    dzdξdκ = -im / (1 - im * ϵ) * ξ

    return (hypgeom_u_dzda(a, b, z) * dadκ + hypgeom_u_dz(a, b, z, 2) * dzdκ) * dzdξ +
           hypgeom_u_dz(a, b, z) * dzdξdκ
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
    dzdξ = -im * κ / (1 - im * ϵ) * ξ

    return exp(z) * (hypgeom_u(b - a, b, -z) - hypgeom_u_dz(b - a, b, -z)) * dzdξ
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
