"""
    solution_zero(μ, κ, ϵ, ξ₁, λ::CGLParams)

Let `Q` be the solution to [`ivp_zero_complex`](@ref). This function
computes `[Q(ξ₁), d(Q)(ξ₁)]`.
"""
function solution_zero(μ::Arb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    sol = solution_zero_capd(μ, κ, ϵ, ξ₁, λ)
    return SVector(Acb(sol[1], sol[2]), Acb(sol[3], sol[4]))
end

function solution_zero(μ, κ, ϵ, ξ₁, λ::CGLParams)
    sol = solution_zero_float(μ, κ, ϵ, ξ₁, λ)
    return SVector(complex(sol[1], sol[2]), complex(sol[3], sol[4]))
end

"""
    solution_zero_jacobian_kappa(μ, κ, ϵ, ξ₁, λ::CGLParams)

Let `Q` be the solution to [`ivp_zero_complex`](@ref). This function
computes `[Q(ξ₁), d(Q)(ξ₁)]` as well as the Jacobian w.r.t. the
parameters `μ` and `κ`. The Jacobian is given by
```
[
d(Q(ξ₁), μ) d(Q(ξ₁), κ)
d(d(Q)(ξ₁), μ) d((Q)(ξ₁), κ)
]
```
where we use `d(Q, μ)` to denote the derivative of `Q` w.r.t. `μ`.
"""
function solution_zero_jacobian_kappa(μ::Arb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    J = solution_zero_jacobian_kappa_capd(μ, κ, ϵ, ξ₁, λ)

    return SMatrix{2,2}(
        Acb(J[1, 1], J[2, 1]),
        Acb(J[3, 1], J[4, 1]),
        Acb(J[1, 2], J[2, 2]),
        Acb(J[3, 2], J[4, 2]),
    )
end

function solution_zero_jacobian_kappa(
    μ::Float64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    J = solution_zero_jacobian_kappa_float(μ, κ, ϵ, ξ₁, λ)

    return SMatrix{2,2}(
        complex(J[1, 1], J[2, 1]),
        complex(J[3, 1], J[4, 1]),
        complex(J[1, 2], J[2, 2]),
        complex(J[3, 2], J[4, 2]),
    )
end

"""
    solution_zero_jacobian_epsilon(μ, κ, ξ₁, λ::CGLParams)

Let `Q` be the solution to [`ivp_zero_complex`](@ref). This function
computes `[Q(ξ₁), d(Q)(ξ₁)]` as well as the Jacobian w.r.t. the
parameters `μ` and `ϵ`. The Jacobian is given by
```
[
d(Q(ξ₁), μ) d(Q(ξ₁), ϵ)
d(d(Q)(ξ₁), μ) d((Q)(ξ₁), ϵ)
]
```
where we use `d(Q, μ)` to denote the derivative of `Q` w.r.t. `μ`.
"""
function solution_zero_jacobian_epsilon(μ::Arb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    J = solution_zero_jacobian_epsilon_capd(μ, κ, ϵ, ξ₁, λ)

    return SMatrix{2,2}(
        Acb(J[1, 1], J[2, 1]),
        Acb(J[3, 1], J[4, 1]),
        Acb(J[1, 2], J[2, 2]),
        Acb(J[3, 2], J[4, 2]),
    )
end

function solution_zero_jacobian_epsilon(
    μ::Float64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    J = solution_zero_jacobian_epsilon_float(μ, κ, ϵ, ξ₁, λ)

    return SMatrix{2,2}(
        complex(J[1, 1], J[2, 1]),
        complex(J[3, 1], J[4, 1]),
        complex(J[1, 2], J[2, 2]),
        complex(J[3, 2], J[4, 2]),
    )
end
