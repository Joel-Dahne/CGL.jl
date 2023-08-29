"""
    solution_zero(γ, κ, ξ₁, λ::AbstractGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes `[Q(ξ₁), d(Q)(ξ₁)]`.
"""
function solution_infinity(γ::Acb, κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    # FIXME
    I_E_bound = Acb(0)
    I_E_dξ_bound = Acb(0)

    # FIXME
    I_P_bound = Acb(0)
    I_P_dξ_bound = Acb(0)

    Q = P(ξ₁, (λ, κ)) * (γ + I_E_bound) + E(ξ₁, (λ, κ)) * I_P_bound

    dQ =
        P(ξ₁, (λ, κ)) * I_E_dξ_bound +
        P_dξ(ξ₁, (λ, κ)) * (γ + I_E_bound) +
        E(ξ₁, (λ, κ)) * I_P_dξ_bound +
        E_dξ(ξ₁, (λ, κ)) * I_P_bound

    return SVector(Q, dQ)
end

function solution_infinity(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::AbstractGLParams{Float64},
)
    # FIXME: Take into account the rest of Q for this case as well?
    Q = γ * P(ξ₁, (λ, κ))
    dQ = γ * P_dξ(ξ₁, (λ, κ))

    return SVector(Q, dQ)
end


"""
    solution_infinity_jacobian(γ, κ, ξ₁, λ::AbstractGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes Jacobian w.r.t. the parameters `γ` and `κ` of
`[Q(ξ₁), d(Q)(ξ₁)]. The Jacobian is given by
```
[
d(Q(ξ₁), μ) d(Q(ξ₁), κ)
d(d(Q)(ξ₁), μ) d((Q)(ξ₁), κ)
]
```
where we use `d(Q, μ)` to denote the derivative of `Q` w.r.t. `μ`.
"""
function solution_infinity_jacobian(γ::Acb, κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    # FIXME
    I_E_bound = Acb(0)
    I_E_dξ_bound = Acb(0)
    I_E_dγ_bound = Acb(0)
    I_E_dκ_bound = Acb(0)
    I_E_dξ_dγ_bound = Acb(0)
    I_E_dξ_dκ_bound = Acb(0)

    # FIXME
    I_P_bound = Acb(0)
    I_P_dξ_bound = Acb(0)
    I_P_dγ_bound = Acb(0)
    I_P_dκ_bound = Acb(0)
    I_P_dξ_dγ_bound = Acb(0)
    I_P_dξ_dκ_bound = Acb(0)

    Q_dγ = P(ξ₁, (λ, κ)) * (one(γ) + I_E_dγ_bound) + E(ξ₁, (λ, κ)) * I_P_dγ_bound

    dQ_dγ =
        P(ξ₁, (λ, κ)) * I_E_dξ_dγ_bound +
        P_dξ(ξ₁, (λ, κ)) * (one(γ) + I_E_dγ_bound) +
        E(ξ₁, (λ, κ)) * I_P_dξ_dγ_bound +
        E_dξ(ξ₁, (λ, κ)) * I_P_dγ_bound

    Q_dκ =
        P(ξ₁, (λ, κ)) * I_E_dκ_bound +
        P_dκ(ξ₁, (λ, κ)) * (γ + I_E_bound) +
        E(ξ₁, (λ, κ)) * I_P_dκ_bound +
        E_dκ(ξ₁, (λ, κ)) * I_P_bound

    dQ_dκ =
        P(ξ₁, (λ, κ)) * I_E_dξ_dκ_bound +
        P_dξ(ξ₁, (λ, κ)) * I_E_dκ_bound +
        P_dκ(ξ₁, (λ, κ)) * I_E_dξ_bound +
        P_dξ_dκ(ξ₁, (λ, κ)) * (γ + I_E_bound) +
        E(ξ₁, (λ, κ)) * I_P_dξ_dκ_bound +
        E_dξ(ξ₁, (λ, κ)) * I_P_dκ_bound +
        E_dκ(ξ₁, (λ, κ)) * I_P_dξ_bound +
        E_dξ_dκ(ξ₁, (λ, κ)) * I_P_bound

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end

function solution_infinity_jacobian(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::AbstractGLParams{Float64},
)
    # FIXME: Take into account the rest of Q
    Q_dγ = complex(P(ξ₁, (λ, κ)))
    dQ_dγ = complex(P_dξ(ξ₁, (λ, κ)))

    Q_dκ = γ * P_dκ(ξ₁, (λ, κ))
    dQ_dκ = γ * P_dξ_dκ(ξ₁, (λ, κ))

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end
