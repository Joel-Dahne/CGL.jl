export GLParams

"""
    GLParams{T}(d, ω, σ, ϵ, δ)
"""
struct GLParams{T}
    d::Int
    ω::T
    σ::T
    ϵ::T
    δ::T
end
