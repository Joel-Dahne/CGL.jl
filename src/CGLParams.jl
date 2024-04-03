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

function sverak_params(
    T::Type{Float64},
    i::Integer = 1,
    d::Integer = 1;
    ξ₁::Union{Float64,Nothing} = nothing,
)
    # Initial approximation from https://doi.org/10.1002/cpa.3006
    if d == 1
        ϵ = 0.0
        λ = CGLParams{T}(1, 1.0, 2.3, 0.0)

        μs = [1.23204, 0.78308, 1.12389, 0.88393, 1.07969, 0.92761, 1.05707, 0.94914]
        κs = [0.85310, 0.49322, 0.34680, 0.26678, 0.21621, 0.18192, 0.15667, 0.13749]
        ξ₁s = Float64[20, 20, 20, 25, 40, 50, 55, 90]
    elseif d == 3
        ϵ = 0.0
        λ = CGLParams{T}(3, 1.0, 1.0, 0.0)

        μs = [1.88619, 0.84142, 1.10919, 0.94337, 1.01123]
        κs = [0.91710, 0.32129, 0.22259, 0.16961, 0.13738]
        ξ₁s = Float64[60, 60, 60, 60, 60]
    else
        error("only contains values d = 1 or d = 3")
    end

    ξ₁ = something(ξ₁, ξ₁s[i])

    # Refine the approximation
    μ, γ, κ = refine_approximation(μs[i], κs[i], ϵ, ξ₁, λ)

    return μ, γ, κ, ϵ, ξ₁, λ
end

function sverak_params(T::Type, i::Integer = 1, d::Integer = 1; ξ₁ = nothing)
    μ, γ, κ, ϵ, ξ₁, λ =
        sverak_params(Float64, i, d; ξ₁ = convert(Union{Float64,Nothing}, ξ₁))

    if T == Arb
        return T(μ), Acb(γ), T(κ), T(ϵ), T(ξ₁), CGLParams{T}(λ)
    else
        return T(μ), complex(T)(γ), T(κ), T(ϵ), T(ξ₁), CGLParams{T}(λ)
    end
end
