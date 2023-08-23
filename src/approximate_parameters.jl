"""
    approximate_parameters_F(ξ₁::Float64, λ::AbstractGLParams{Float64})

Return a function `F(μ, κ)` for computing the real and imaginary part
of
```
d(Q)(ξ₁) - γ * P_dξ(ξ₁, (p, κ))
```
where
```
Q, d(Q) = solution_zero(μ, κ, ξ₁, λ)
γ = Q(ξ₁) / P(ξ₁, (p, κ))
```
"""
function approximate_parameters_F(ξ₁::Float64, λ::AbstractGLParams{Float64})
    return (μ, κ) -> begin
        Q, dQ = solution_zero(μ, κ, ξ₁, λ)

        γ = Q / P_asym_approx(ξ₁, (λ, κ))

        res = dQ - γ * P_dξ_asym_approx(ξ₁, (λ, κ))

        return SVector(real(res), imag(res))
    end
end

"""
    approximate_parameters(μ₀::T, κ₀::T, ξ₁::T, λ::AbstractGLParams{T}) where {T}

Compute `μ, γ, κ` such that
```
d(Q)(ξ₁) = γ * P_dξ(ξ₁, (p, κ))
```
where
```
Q, d(Q) = solution_zero(μ, κ, ξ₁, λ)
γ = Q(ξ₁) / P(ξ₁, (p, κ))
```
"""
function approximate_parameters(μ₀::T, κ₀::T, ξ₁::T, λ::AbstractGLParams{T}) where {T}
    F = splat(approximate_parameters_F(ξ₁, λ))
    sol = nlsolve(F, [μ₀, κ₀], xtol = 1e-12, ftol = 1e-12)
    μ, κ = sol.zero

    γ = let (Q, dQ) = solution_zero(μ, κ, ξ₁, λ)
        Q / P(ξ₁, (λ, κ))
    end

    return μ, γ, κ
end
