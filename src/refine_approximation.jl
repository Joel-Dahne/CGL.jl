"""
    refine_approximation(μ₀, γ₀, κ₀, ϵ, ξ₁, λ; return_convergence, verbose)
    refine_approximation(μ₀, κ₀, ϵ, ξ₁, λ; return_convergence, verbose)

Given an initial approximation to a zero of
```
G(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)
```
compute a refined approximation. The version without `γ₀` uses the
initial approximation
```
γ₀ = solution_zero(μ₀, κ₀, ϵ₀, ξ₁, λ)[1] / P(ξ₁, κ₀, ϵ, λ)
```
for it.

If the argument `return_convergence` can bet set to `Val{true}()` then
the first value returned is a flag indicating if the refinement seems
to have converged or not.
"""
function refine_approximation(
    μ₀::Float64,
    γ₀::ComplexF64,
    κ₀::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64};
    return_convergence::Union{Val{false},Val{true}} = Val{false}(),
    verbose = false,
)
    F = x -> G(x..., ϵ, ξ₁, λ)
    # 100 iterations should be enough to saturate the convergence in
    # practice
    sol = nlsolve(F, [μ₀, real(γ₀), imag(γ₀), κ₀], iterations = 100, ftol = 1e-10)

    converged = sol.residual_norm < 1e-6

    if verbose && !converged
        @warn "Low precision when refining approximation" ϵ sol.residual_norm
    end

    μ, γ_real, γ_imag, κ = sol.zero

    if return_convergence isa Val{false}
        return μ, complex(γ_real, γ_imag), κ
    else
        return converged, μ, complex(γ_real, γ_imag), κ
    end
end

function refine_approximation(
    μ₀::Arb,
    γ₀::Acb,
    κ₀::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    return_convergence::Union{Val{false},Val{true}} = Val{false}(),
    extra_newton = 2,
    verbose = false,
)
    converged, μ, γ, κ = refine_approximation(
        Float64(μ₀),
        ComplexF64(γ₀),
        Float64(κ₀),
        Float64(ϵ),
        Float64(ξ₁),
        CGLParams{Float64}(λ),
        return_convergence = Val{true}();
        verbose,
    )

    x = SVector{4,Arb}(μ, real(γ), imag(γ), κ)

    # Potentially do a few Newton iterations with Arb.
    for _ = 1:(converged*extra_newton)
        y = G(x..., midpoint(Arb, ϵ), ξ₁, λ)

        # If the enclosure already contains zero more iterations are
        # unlikely to help much.
        all(Arblib.contains_zero, y) && break

        J = G_jacobian_kappa(x..., midpoint(Arb, ϵ), ξ₁, λ)

        # We don't need high precision here, so it is fine to do the
        # update in Float64. In rare cases J is (numerically) singular
        # and we catch the exception then.
        x = try
            Arb.(Float64.(x) - Float64.(J) \ Float64.(y))
        catch e
            e isa SingularException || rethrow(e)
            x
        end
    end

    if return_convergence isa Val{false}
        return x[1], Acb(x[2], x[3]), x[4]
    else
        return converged, x[1], Acb(x[2], x[3]), x[4]
    end
end

function refine_approximation(
    μ₀::Float64,
    κ₀::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64};
    return_convergence::Union{Val{false},Val{true}} = Val{false}(),
    verbose = false,
)
    γ₀ = solution_zero(μ₀, κ₀, ϵ, ξ₁, λ)[1] / P(ξ₁, κ₀, ϵ, λ)

    return refine_approximation(μ₀, γ₀, κ₀, ϵ, ξ₁, λ; return_convergence, verbose)
end

function refine_approximation(
    μ₀::Arb,
    κ₀::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    return_convergence::Union{Val{false},Val{true}} = Val{false}(),
    extra_newton = 2,
    verbose = false,
)
    μ₀_F64 = Float64(μ₀)
    κ₀_F64 = Float64(κ₀)
    ϵ_F64 = Float64(ϵ)
    ξ₁_F64 = Float64(ξ₁)
    λ_F64 = CGLParams{Float64}(λ)

    γ₀_F64 =
        solution_zero(μ₀_F64, κ₀_F64, ϵ_F64, ξ₁_F64, λ_F64)[1] /
        P(ξ₁_F64, κ₀_F64, ϵ_F64, λ_F64)

    return refine_approximation(
        μ₀,
        Acb(γ₀_F64),
        κ₀,
        ϵ,
        ξ₁,
        λ;
        return_convergence,
        extra_newton,
        verbose,
    )
end

"""
    refine_approximation_with_interpolation((μ₁, μ₂), (κ₁, κ₂), (ϵ₁, ϵ₂), ϵ, ξ₁, λ)

Compute a refined approximation to a zero of
```
G(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)
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
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    verbose = false,
) where {T}
    t = if T == Arb
        @assert ϵ₁ <= midpoint(ϵ) <= ϵ₂ # Not needed but in practice the case
        (ϵ₂ - midpoint(ϵ)) / (ϵ₂ - ϵ₁)
    else
        @assert ϵ₁ <= ϵ <= ϵ₂  # Not needed but in practice the case
        (ϵ₂ - ϵ) / (ϵ₂ - ϵ₁)
    end

    μ₀ = (1 - t) * μ₁ + t * μ₂
    κ₀ = (1 - t) * κ₁ + t * κ₂

    return refine_approximation(μ₀, κ₀, ϵ, ξ₁, λ)
end
