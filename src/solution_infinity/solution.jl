"""
    solution_zero(μ, κ, ξ₁, λ::AbstractGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes `[Q(ξ₁), d(Q)(ξ₁)]`.
"""
function solution_infinity(γ::Acb, κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    # FIXME: Take into account the rest of Q
    Q = γ * P(ξ₁, (λ, κ))
    dQ = γ * P_dξ(ξ₁, (λ, κ))

    return SVector(Q, dQ)
end

function solution_infinity(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::AbstractGLParams{Float64},
)
    # FIXME: Take into account the rest of Q
    Q = γ * P(ξ₁, (λ, κ))
    dQ = γ * P_dξ(ξ₁, (λ, κ))

    return SVector(Q, dQ)
end

"""
    solution_infinity_jacobian(μ, κ, ξ₁, λ::AbstractGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes `[Q(ξ₁), d(Q)(ξ₁)]` as well as the Jacobian w.r.t.
the parameters `γ` and `κ`. The Jacobian is given by
```
[
d(Q(ξ₁), μ) d(Q(ξ₁), κ)
d(d(Q)(ξ₁), μ) d((Q)(ξ₁), κ)
]
```
where we use `d(Q, μ)` to denote the derivative of `Q` w.r.t. `μ`.
"""
function solution_infinity_jacobian(γ::Acb, κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    # FIXME: Take into account the rest of Q
    Q = γ * P(ξ₁, (λ, κ))
    dQ = γ * P_dξ(ξ₁, (λ, κ))

    Q_dγ = Acb(P(ξ₁, (λ, κ)))
    dQ_dγ = Acb(P_dξ(ξ₁, (λ, κ)))

    Q_dκ = γ * P_dκ(ξ₁, (λ, κ))
    dQ_dκ = γ * P_dξdκ(ξ₁, (λ, κ))

    jacobian = SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)

    return SVector(Q, dQ), jacobian
end

function solution_infinity_jacobian(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::AbstractGLParams{Float64},
)
    # FIXME: Take into account the rest of Q
    Q = γ * P(ξ₁, (λ, κ))
    dQ = γ * P_dξ(ξ₁, (λ, κ))

    Q_dγ = complex(P(ξ₁, (λ, κ)))
    dQ_dγ = complex(P_dξ(ξ₁, (λ, κ)))

    Q_dκ = γ * P_dκ(ξ₁, (λ, κ))
    dQ_dκ = γ * P_dξdκ(ξ₁, (λ, κ))

    jacobian = SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)

    return SVector(Q, dQ), jacobian
end
