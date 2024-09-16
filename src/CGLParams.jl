export CGLParams

"""
    CGLParams{T}(d, ω, σ, δ)
"""
struct CGLParams{T}
    d::Int
    ω::T
    σ::T
    δ::T
end

CGLParams{T}(λ::CGLParams; d = λ.d, ω = λ.ω, σ = λ.σ, δ = λ.δ) where {T} =
    CGLParams{T}(d, ω, σ, δ)

CGLParams(λ::CGLParams{T}; d = λ.d, ω = λ.ω, σ = λ.σ, δ = λ.δ) where {T} =
    CGLParams{T}(d, ω, σ, δ)

Base.isequal(λ1::CGLParams{T}, λ2::CGLParams{T}) where {T} =
    isequal(λ1.d, λ2.d) && isequal(λ1.ω, λ2.ω) && isequal(λ1.σ, λ2.σ) && isequal(λ1.δ, λ2.δ)

"""
    scale_params(μ, γ, κ, ϵ, ξ₁, λ::CGLParams; scaling)
    scale_params(μ, γ_real, γ_imag, κ, ϵ, ξ₁, λ::CGLParams; scaling)
    scale_params(μ, γ, κ, λ; scaling)
    scale_params(μ, γ_real, γ_imag, κ, λ; scaling)
    scale_params(μ_γ_κ::SVector{4}, λ; scaling)

The equation has a scaling symmetry for `scaling > 0`. This return the
scaled parameters.

The version taking only `μ`, `γ` and `κ` is useful for scaling back
the output of [`G_solve`](@ref).
"""
function scale_params(μ, γ, κ, ϵ, ξ₁, λ::CGLParams; scaling)
    μ_scaled, γ_scaled, κ_scaled = scale_params(μ, γ, κ, λ; scaling)
    ϵ_scaled = ϵ # No scaling for ϵ
    ξ₁_scaled = ξ₁ / scaling
    λ_scaled = CGLParams(λ, ω = λ.ω * scaling^2)

    return μ_scaled, γ_scaled, κ_scaled, ϵ_scaled, ξ₁_scaled, λ_scaled
end

function scale_params(μ, γ_real, γ_imag, κ, ϵ, ξ₁, λ::CGLParams; scaling)
    μ_scaled, γ_real_scaled, γ_imag_scaled, κ_scaled =
        scale_params(μ, γ_real, γ_imag, κ, λ; scaling)
    ϵ_scaled = ϵ # No scaling for ϵ
    ξ₁_scaled = ξ₁ / scaling
    λ_scaled = CGLParams(λ, ω = λ.ω * scaling^2)

    return μ_scaled, γ_real_scaled, γ_imag_scaled, κ_scaled, ϵ_scaled, ξ₁_scaled, λ_scaled
end

function scale_params(μ, γ, κ, λ::CGLParams; scaling)
    μ_scaled = μ * scaling^(1 / λ.σ)
    γ_scaled = γ * scaling^(1 / λ.σ)
    κ_scaled = κ * scaling^2

    return μ_scaled, γ_scaled, κ_scaled
end

function scale_params(μ, γ_real, γ_imag, κ, λ::CGLParams; scaling)
    μ_scaled = μ * scaling^(1 / λ.σ)
    γ_real_scaled = γ_real * scaling^(1 / λ.σ)
    γ_imag_scaled = γ_imag * scaling^(1 / λ.σ)
    κ_scaled = κ * scaling^2

    return μ_scaled, γ_real_scaled, γ_imag_scaled, κ_scaled
end

scale_params(μ_γ_κ::SVector{4}, λ::CGLParams; scaling) =
    SVector{4}(scale_params(μ_γ_κ..., λ; scaling))

function sverak_params(
    ::Type{T},
    j::Integer = 1,
    d::Integer = 1;
    ξ₁::Union{Real,Nothing} = nothing,
) where {T}
    # Initial approximation from https://doi.org/10.1002/cpa.3006
    if d == 1
        ϵ = T(0.0)
        λ = CGLParams{T}(1, 1.0, 2.3, 0.0)

        μs = T[1.23204, 0.78308, 1.12389, 0.88393, 1.07969, 0.92761, 1.05707, 0.94914]
        κs = T[0.85310, 0.49322, 0.34680, 0.26678, 0.21621, 0.18192, 0.15667, 0.13749]
        ξ₁s = T[20, 20, 20, 30, 40, 50, 60, 90]
    elseif d == 3
        ϵ = T(0.0)
        λ = CGLParams{T}(3, 1.0, 1.0, 0.0)

        μs = T[1.88565, 0.84142, 1.10919, 0.94337, 1.01123]
        κs = T[0.91734, 0.32129, 0.22259, 0.16961, 0.13738]
        ξ₁s = T[15, 25, 25, 40, 40]
    else
        error("only contains values d = 1 or d = 3")
    end

    ξ₁ = convert(T, something(ξ₁, ξ₁s[j]))

    # Refine the approximation
    μ, γ, κ = refine_approximation_fix_epsilon(μs[j], κs[j], ϵ, ξ₁, λ)

    return μ, γ, κ, ϵ, ξ₁, λ
end
