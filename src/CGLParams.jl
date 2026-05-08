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

CGLParams{T}(Λ::CGLParams; d = Λ.d, ω = Λ.ω, σ = Λ.σ, δ = Λ.δ) where {T} =
    CGLParams{T}(d, ω, σ, δ)

CGLParams(Λ::CGLParams{T}; d = Λ.d, ω = Λ.ω, σ = Λ.σ, δ = Λ.δ) where {T} =
    CGLParams{T}(d, ω, σ, δ)

Base.isequal(Λ1::CGLParams{T}, Λ2::CGLParams{T}) where {T} =
    isequal(Λ1.d, Λ2.d) && isequal(Λ1.ω, Λ2.ω) && isequal(Λ1.σ, Λ2.σ) && isequal(Λ1.δ, Λ2.δ)

"""
    scale_params(μ, γ, κ, ϵ, ξ₁, Λ::CGLParams; scaling)
    scale_params(μ, γ_real, γ_imag, κ, ϵ, ξ₁, Λ::CGLParams; scaling)
    scale_params(μ, γ, κ, Λ; scaling)
    scale_params(μ, γ_real, γ_imag, κ, Λ; scaling)
    scale_params(μ_γ_κ::SVector{4}, Λ; scaling)

The equation has a scaling symmetry for `scaling > 0`. This return the
scaled parameters.

The version taking only `μ`, `γ` and `κ` is useful for scaling back
the output of [`G_solve`](@ref).
"""
function scale_params(μ, γ, κ, ϵ, ξ₁, Λ::CGLParams; scaling)
    μ_scaled, γ_scaled, κ_scaled = scale_params(μ, γ, κ, Λ; scaling)
    ϵ_scaled = ϵ # No scaling for ϵ
    ξ₁_scaled = ξ₁ / scaling
    Λ_scaled = CGLParams(Λ, ω = Λ.ω * scaling^2)

    return μ_scaled, γ_scaled, κ_scaled, ϵ_scaled, ξ₁_scaled, Λ_scaled
end

function scale_params(μ, γ_real, γ_imag, κ, ϵ, ξ₁, Λ::CGLParams; scaling)
    μ_scaled, γ_real_scaled, γ_imag_scaled, κ_scaled =
        scale_params(μ, γ_real, γ_imag, κ, Λ; scaling)
    ϵ_scaled = ϵ # No scaling for ϵ
    ξ₁_scaled = ξ₁ / scaling
    Λ_scaled = CGLParams(Λ, ω = Λ.ω * scaling^2)

    return μ_scaled, γ_real_scaled, γ_imag_scaled, κ_scaled, ϵ_scaled, ξ₁_scaled, Λ_scaled
end

function scale_params(μ, γ, κ, Λ::CGLParams; scaling)
    μ_scaled = μ * scaling^(1 / Λ.σ)
    γ_scaled = γ * scaling^(1 / Λ.σ)
    κ_scaled = κ * scaling^2

    return μ_scaled, γ_scaled, κ_scaled
end

function scale_params(μ, γ_real, γ_imag, κ, Λ::CGLParams; scaling)
    μ_scaled = μ * scaling^(1 / Λ.σ)
    γ_real_scaled = γ_real * scaling^(1 / Λ.σ)
    γ_imag_scaled = γ_imag * scaling^(1 / Λ.σ)
    κ_scaled = κ * scaling^2

    return μ_scaled, γ_real_scaled, γ_imag_scaled, κ_scaled
end

scale_params(μ_γ_κ::SVector{4}, Λ::CGLParams; scaling) =
    SVector{4}(scale_params(μ_γ_κ..., Λ; scaling))

function sverak_params(
    ::Type{T},
    j::Integer = 1,
    d::Integer = 1;
    ξ₁::Union{Real,Nothing} = nothing,
    ξ₁_for_branch::Bool = false,
) where {T}
    # Initial approximation from https://doi.org/10.1002/cpa.3006
    if d == 1
        ϵ = T(0.0)

        # Make sure we take an enclosure of 2.3 in the case of T =
        # Arb.
        σ = T == Arb ? Arb("2.3") : T(2.3)

        Λ = CGLParams{T}(1, 1.0, σ, 0.0)

        μs = T[1.23204, 0.78308, 1.12389, 0.88393, 1.07969, 0.92761, 1.05707, 0.94914]
        κs = T[0.85310, 0.49322, 0.34680, 0.26678, 0.21621, 0.18192, 0.15667, 0.13749]
        if ξ₁_for_branch
            ξ₁s = T[20, 20, 20, 30, 40, 50, 60, 90]
        else
            ξ₁s = T[10, 15, 20, 25, 40, 45, 60, 75]
        end
    elseif d == 3
        ϵ = T(0.0)
        Λ = CGLParams{T}(3, 1.0, 1.0, 0.0)

        μs = T[1.88565, 0.84142, 1.10919, 0.94337, 1.01123]
        κs = T[0.91734, 0.32129, 0.22259, 0.16961, 0.13738]
        if ξ₁_for_branch
            ξ₁s = T[15, 25, 25, 40, 40]
        else
            ξ₁s = T[60, 140, 60, 60, 60]
        end
    else
        error("only contains values d = 1 or d = 3")
    end

    ξ₁ = convert(T, something(ξ₁, ξ₁s[j]))

    # Refine the approximation
    μ, γ, κ = refine_approximation_fix_epsilon(μs[j], κs[j], ϵ, ξ₁, Λ)

    return μ, γ, κ, ϵ, ξ₁, Λ
end
