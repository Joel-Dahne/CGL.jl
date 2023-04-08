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
