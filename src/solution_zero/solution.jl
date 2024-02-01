"""
    solution_zero(μ, κ, ξ₁, λ::CGLParams)

Let `Q` be the solution to [`ivp_zero_complex`](@ref). This function
computes `[Q(ξ₁), d(Q)(ξ₁)]`.
"""
function solution_zero(μ::Arb, κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    sol = solution_zero_capd(μ, κ, ξ₁, λ)
    return SVector(Acb(sol[1], sol[2]), Acb(sol[3], sol[4]))
end

function solution_zero(μ, κ, ξ₁, λ::CGLParams)
    sol = solution_zero_float(μ, κ, ξ₁, λ)
    return SVector(complex(sol[1], sol[2]), complex(sol[3], sol[4]))
end

"""
    solution_zero_jacobian(μ, κ, ξ₁, λ::CGLParams)

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
function solution_zero_jacobian(μ::Arb, κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    sol, sol_jac = solution_zero_jacobian_capd(μ, κ, ξ₁, λ)

    res = SVector(Acb(sol[1], sol[2]), Acb(sol[3], sol[4]))

    res_jac = SMatrix{2,2}(
        Acb(sol_jac[1, 1], sol_jac[2, 1]),
        Acb(sol_jac[3, 1], sol_jac[4, 1]),
        Acb(sol_jac[1, 2], sol_jac[2, 2]),
        Acb(sol_jac[3, 2], sol_jac[4, 2]),
    )

    return res, res_jac
end

function solution_zero_jacobian(μ::Float64, κ::Float64, ξ₁::Float64, λ::CGLParams{Float64})
    sol, sol_jac = solution_zero_jacobian_float(μ, κ, ξ₁, λ)

    res = SVector(complex(sol[1], sol[2]), complex(sol[3], sol[4]))

    res_jac = SMatrix{2,2}(
        complex(sol_jac[1, 1], sol_jac[2, 1]),
        complex(sol_jac[3, 1], sol_jac[4, 1]),
        complex(sol_jac[1, 2], sol_jac[2, 2]),
        complex(sol_jac[3, 2], sol_jac[4, 2]),
    )

    return res, res_jac
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
function solution_zero_jacobian_epsilon(μ::Arb, κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    sol, sol_jac = solution_zero_jacobian_epsilon_capd(μ, κ, ξ₁, λ)

    res = SVector(Acb(sol[1], sol[2]), Acb(sol[3], sol[4]))

    res_jac = SMatrix{2,2}(
        Acb(sol_jac[1, 1], sol_jac[2, 1]),
        Acb(sol_jac[3, 1], sol_jac[4, 1]),
        Acb(sol_jac[1, 2], sol_jac[2, 2]),
        Acb(sol_jac[3, 2], sol_jac[4, 2]),
    )

    return res, res_jac
end
