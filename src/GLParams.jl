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
    ω, σ, ϵ, δ = promote(ω, σ, ϵ, δ)
    if isone(d)
        return GLParamsd1{typeof(ω)}(ω, σ, ϵ, δ)
    else
        return GLParams{typeof(ω)}(d, ω, σ, ϵ, δ)
    end
end

Base.convert(S::Type{<:AbstractGLParams{T}}, p::AbstractGLParams) where {T} =
    S(p.d, p.ω, p.σ, p.ϵ, p.δ)
