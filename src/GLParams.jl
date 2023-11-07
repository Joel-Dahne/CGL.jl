export AbstractGLParams, GLParams, GLParamsd1, gl_params

abstract type AbstractGLParams{T} end

"""
    GLParams{T}(d, ω, σ, ϵ, δ)
"""
struct GLParams{T} <: AbstractGLParams{T}
    d::Int
    ω::T
    σ::T
    ϵ::T
    δ::T
end

"""
    GLParamsd1{T}(ω, σ, ϵ, δ)

Like [`GLParams`](@ref) but with `d` fixed to `1`.
"""
struct GLParamsd1{T} <: AbstractGLParams{T}
    ω::T
    σ::T
    ϵ::T
    δ::T
end

function GLParamsd1{T}(d, ω, σ, ϵ, δ) where {T}
    isone(d) || throw(ArgumentError("only accepts d == 1"))
    return GLParamsd1{T}(ω, σ, ϵ, δ)
end

function Base.getproperty(p::GLParamsd1, name::Symbol)
    if name == :d
        return 1
    else
        return getfield(p, name)
    end
end

function gl_params(d::Int, ω, σ, ϵ, δ)
    T = promote_type(typeof(ω), typeof(σ), typeof(ϵ), typeof(δ))
    return gl_params(T, d, ω, σ, ϵ, δ)
end

function gl_params(T::Type, d::Int, ω, σ, ϵ, δ)
    ω, σ, ϵ, δ = convert(T, ω), convert(T, σ), convert(T, ϵ), convert(T, δ)
    if isone(d)
        return GLParamsd1{T}(ω, σ, ϵ, δ)
    else
        return GLParams{T}(d, ω, σ, ϵ, δ)
    end
end

gl_params(p::AbstractGLParams) = gl_params(p.d, p.ω, p.σ, p.ϵ, p.δ)
gl_params(T::Type, p::AbstractGLParams) = gl_params(T, p.d, p.ω, p.σ, p.ϵ, p.δ)

function _params(T::Type{Float64}, i::Integer = 1, d::Integer = 1; ξ₁::Float64 = 30.0)
    if d == 1
        λ = gl_params(1, 1.0, 2.3, 0.0, 0.0)

        μs = [1.23204, 0.78308, 1.12389, 0.88393, 1.07969, 0.92761, 1.05707, 0.94914]
        κs = [0.85310, 0.49322, 0.34680, 0.26678, 0.21621, 0.18192, 0.15667, 0.13749]

    elseif d == 3
        λ = gl_params(T, 3, 1.0, 1.0, 0.0, 0.0)

        μs = [1.88619, 0.84142, 1.10919, 0.94337, 1.01123]
        κs = [0.91710, 0.32129, 0.22259, 0.16961, 0.13738]
    else
        error("only contains values d = 1 or d = 3")
    end

    # Compute a first approximation, giving γ
    μ₀, γ₀, κ₀ = approximate_parameters_simple(μs[i], κs[i], ξ₁, λ)

    # Compute a better approximation
    μ, γ, κ = approximate_parameters(μ₀, γ₀, κ₀, ξ₁, λ)

    return μ, γ, κ, ξ₁, λ
end

function _params(T::Type, i::Integer = 1, d::Integer = 1; ξ₁ = 30.0)
    μ, γ, κ, ξ₁, λ = _params(Float64, i, d; ξ₁ = convert(Float64, ξ₁))

    if T == Arb
        return T(μ), Acb(γ), T(κ), T(ξ₁), gl_params(T, λ)
    else
        return T(μ), complex(T)(γ), T(κ), T(ξ₁), gl_params(T, λ)
    end
end
