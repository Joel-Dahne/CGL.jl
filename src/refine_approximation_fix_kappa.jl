function refine_approximation_fix_kappa(
    μ₀::Float64,
    γ₀::ComplexF64,
    κ::Float64,
    ϵ₀::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64};
    verbose = false,
)
    F = x -> G(x[1:3]..., κ, x[4], ξ₁, λ)
    # 100 iterations should be enough to saturate the convergence in
    # practice
    sol = nlsolve(F, [μ₀, real(γ₀), imag(γ₀), ϵ₀], iterations = 100, ftol = 1e-10)

    if verbose && sol.residual_norm > 1e-8
        @warn "Low precision when refining approximation" λ.ϵ sol.residual_norm
    end

    μ, γ_real, γ_imag, ϵ = sol.zero

    return μ, complex(γ_real, γ_imag), ϵ
end

function refine_approximation_fix_kappa(
    μ₀::Arb,
    γ₀::Acb,
    κ::Arb,
    ϵ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    extra_newton = 2,
    verbose = false,
)
    μ, γ, ϵ = refine_approximation_fix_kappa(
        Float64(μ₀),
        ComplexF64(γ₀),
        Float64(κ),
        Float64(ϵ₀),
        Float64(ξ₁),
        CGLParams{Float64}(λ);
        verbose,
    )

    x = SVector{4,Arb}(μ, real(γ), imag(γ), ϵ₀)

    # Potentially do a few Newton iterations with Arb.
    for _ = 1:extra_newton
        y = G(x[1:3]..., midpoint(Arb, κ), x[4], ξ₁, λ)

        # If the enclosure already contains zero more iterations are
        # unlikely to help much.
        all(Arblib.contains_zero, y) && break

        J = G_jacobian_epsilon(x[1:3]..., midpoint(Arb, κ), x[4], ξ₁, λ)

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

    return x[1], Acb(x[2], x[3]), x[4]
end

function refine_approximation_fix_kappa(
    μ₀::Float64,
    κ::Float64,
    ϵ₀::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64};
    verbose = false,
)
    γ₀ = solution_zero(μ₀, κ, ϵ₀, ξ₁, λ)[1] / P(ξ₁, κ, ϵ₀, λ)

    return refine_approximation_fix_kappa(μ₀, γ₀, κ, ϵ₀, ξ₁, λ; verbose)
end

function refine_approximation_fix_kappa(
    μ₀::Arb,
    κ::Arb,
    ϵ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    extra_newton = 2,
    verbose = false,
)
    μ₀_F64 = Float64(μ₀)
    κ_F64 = Float64(κ)
    ϵ₀_F64 = Float64(ϵ₀)
    ξ₁_F64 = Float64(ξ₁)
    λ_F64 = CGLParams{Float64}(λ)

    γ₀_F64 =
        solution_zero(μ₀_F64, κ_F64, ϵ₀_F64, ξ₁_F64, λ_F64)[1] /
        P(ξ₁_F64, κ_F64, ϵ₀_F64, λ_F64)

    return refine_approximation_fix_kappa(
        μ₀,
        Acb(γ₀_F64),
        κ,
        ϵ₀,
        ξ₁,
        λ;
        extra_newton,
        verbose,
    )
end

function refine_approximation_fix_kappa_with_interpolation(
    (μ₁, μ₂)::Tuple{T,T},
    (κ₁, κ₂)::Tuple{T,T},
    (ϵ₁, ϵ₂)::Tuple{T,T},
    κ::T,
    ξ₁::T,
    λ::CGLParams{T};
    verbose = false,
) where {T}
    t = if T == Arb
        @assert κ₂ < midpoint(κ) < κ₁ # Not needed but in practice the case
        (κ₂ - midpoint(κ)) / (κ₂ - κ₁)
    else
        @assert κ₂ < κ < κ₁  # Not needed but in practice the case
        (κ₂ - κ) / (κ₂ - κ₁)
    end

    μ₀ = (1 - t) * μ₁ + t * μ₂
    ϵ₀ = (1 - t) * ϵ₁ + t * ϵ₂

    return refine_approximation_fix_kappa(μ₀, κ, ϵ₀, ξ₁, λ)
end
