"""
    refine_approximation(μ₀, γ₀, κ₀, ξ₁, λ)
    refine_approximation(μ₀, κ₀, ξ₁, λ)

Given an initial approximation to a zero of
```
G_real(μ, real(γ), imag(γ), κ, ξ₁, λ)
```
compute a refined approximation. The version without `γ₀` uses the
initial approximation
```
γ₀ = solution_zero(μ₀, κ₀, ξ₁, λ)[1] / P(ξ₁, (λ, κ₀))
```
for it.
"""
function refine_approximation(
    μ₀::Float64,
    γ₀::ComplexF64,
    κ₀::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    F = x -> G_real(x..., ξ₁, λ)
    # 100 iterations should be enough to saturate the convergence in
    # practice
    sol = nlsolve(F, [μ₀, real(γ₀), imag(γ₀), κ₀], iterations = 100, ftol = 1e-10)
    μ, γ_real, γ_imag, κ = sol.zero

    return μ, complex(γ_real, γ_imag), κ
end

function refine_approximation(μ₀::Arb, γ₀::Acb, κ₀::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    μ, γ, κ = refine_approximation(
        Float64(μ₀),
        ComplexF64(γ₀),
        Float64(κ₀),
        Float64(ξ₁),
        CGLParams{Float64}(λ),
    )

    return Arb(μ), Acb(γ), Arb(κ)
end

function refine_approximation(μ₀::Float64, κ₀::Float64, ξ₁::Float64, λ::CGLParams{Float64})
    γ₀ = solution_zero(μ₀, κ₀, ξ₁, λ)[1] / P(ξ₁, (λ, κ₀))

    return refine_approximation(μ₀, γ₀, κ₀, ξ₁, λ)
end

function refine_approximation(μ₀::Arb, κ₀::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    μ, γ, κ =
        refine_approximation(Float64(μ₀), Float64(κ₀), Float64(ξ₁), CGLParams{Float64}(λ))

    return Arb(μ), Acb(γ), Arb(κ)
end
