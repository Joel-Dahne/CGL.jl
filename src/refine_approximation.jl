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
    λ::CGLParams{Float64};
    verbose = false,
)
    F = x -> G_real(x..., ξ₁, λ)
    # 100 iterations should be enough to saturate the convergence in
    # practice
    sol = nlsolve(F, [μ₀, real(γ₀), imag(γ₀), κ₀], iterations = 100, ftol = 1e-10)

    if verbose && sol.residual_norm > 1e-8
        @warn "Low precision when refining approximation" λ.ϵ sol.residual_norm
    end

    μ, γ_real, γ_imag, κ = sol.zero

    return μ, complex(γ_real, γ_imag), κ
end

function refine_approximation(
    μ₀::Arb,
    γ₀::Acb,
    κ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    verbose = false,
)
    μ, γ, κ = refine_approximation(
        Float64(μ₀),
        ComplexF64(γ₀),
        Float64(κ₀),
        Float64(ξ₁),
        CGLParams{Float64}(λ);
        verbose,
    )

    return Arb(μ), Acb(γ), Arb(κ)
end

function refine_approximation(
    μ₀::Float64,
    κ₀::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64};
    verbose = false,
)
    γ₀ = solution_zero(μ₀, κ₀, ξ₁, λ)[1] / P(ξ₁, (λ, κ₀))

    return refine_approximation(μ₀, γ₀, κ₀, ξ₁, λ; verbose)
end

function refine_approximation(μ₀::Arb, κ₀::Arb, ξ₁::Arb, λ::CGLParams{Arb}; verbose = false)
    μ, γ, κ = refine_approximation(
        Float64(μ₀),
        Float64(κ₀),
        Float64(ξ₁),
        CGLParams{Float64}(λ);
        verbose,
    )

    return Arb(μ), Acb(γ), Arb(κ)
end

"""
    refine_approximation_with_interpolation((μ₁, μ₂), (κ₁, κ₂), (ϵ₁, ϵ₂), ξ₁, λ)

Compute a refined approximation to a zero of
```
G_real(μ, real(γ), imag(γ), κ, ξ₁, λ)
```
The initial approximation for `μ` and `κ` is given by linearly
interpolating between the two input arguments. It is then refined
using [`refine_approximation`](@ref).

The current implementation assumes that `ϵ₁ < λ.ϵ < ϵ₂`, even if this
is not strictly necessary.
"""
function refine_approximation_with_interpolation(
    (μ₁, μ₂)::Tuple{T,T},
    (κ₁, κ₂)::Tuple{T,T},
    (ϵ₁, ϵ₂)::Tuple{T,T},
    ξ₁::T,
    λ::CGLParams{T};
    verbose = false,
) where {T}
    (; ϵ) = λ

    t = if T == Arb
        @assert ϵ₁ < midpoint(ϵ) < ϵ₂ # Not needed but in practice the case
        (ϵ₂ - midpoint(ϵ)) / (ϵ₂ - ϵ₁)
    else
        @assert ϵ₁ < ϵ < ϵ₂  # Not needed but in practice the case
        (ϵ₂ - ϵ) / (ϵ₂ - ϵ₁)
    end

    μ₀ = (1 - t) * μ₁ + t * μ₂
    κ₀ = (1 - t) * κ₁ + t * κ₂

    return refine_approximation(μ₀, κ₀, ξ₁, λ)
end
