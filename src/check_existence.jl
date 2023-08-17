"""
    check_existence_fixed_point(κ::Arb, μ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)

Attempt to prove that the there exists `γ` such that the operator `T`
has a fixed point `u_T(ξ, (p, κ, γ))` satisfying that
```
u_0(ξ₁, (p, κ, μ)) = u_T(ξ₁, (p, κ, γ))
```
Returns `true, γ` on success, where `γ` is a ball enclosing the
existing value of `γ`. On failure it returns false and an
indeterminate value for `γ`.

Here `u_0(ξ, (p, κ, μ))` is the solution satisfying the initial condition
`u_0(0, (p, κ, μ)) = μ`.

# Computing `u_0(ξ₁, (p, κ, μ))`
This is done using [`ode_series_solver`](@ref).

Alternatively if `non_rigorous = true` it uses a non-rigorous solver
from [`DifferentialEquations`](@ref) on the lower and upper bounds of
`κ` and `μ` and takes the union of all four combinations.

# An approximate solution
Using that
```
u_T(ξ, (p, κ, γ)) ≈ γ * P(ξ₁, (p, κ))
```
the equation is approximately satisfied if we take `γ` as
```
u_0(ξ₁, (p, κ, μ)) / P(ξ₁, (p, κ))
```

# Validating the solution
If we let
```
γ₀ = u_0(ξ₁, (p, κ, μ)) / P(ξ₁, (p, κ))
```
be the approximate solution from above we want to prove the existence
of an exact solution nearby this approximation. More precisely we want
to prove the existence of a zero to
```
u_0(ξ₁, (p, κ, μ)) - u_T(ξ₁, (p, κ, γ))
```

We do this using one iteration of the interval Newton method. If we
let `γ_ball` be a ball containing `γ₀` then we need to prove that
```
γ₀ - u_T(ξ₁, (p, κ, γ₀)) / u_T_dγ(ξ₁, (p, κ, γ_ball))
```
is contained in `γ_ball`. Here `u_T_dγ(ξ₁, (p, κ, γ_ball))` denotes
`u_T(ξ₁, (p, κ, γ₀))` differentiated once with respect to `γ`.

## Enclosing `u_T(ξ₁, (p, κ, γ₀))`
If we let
```
C₁, C₂ = C_fix_point(abs(γ₀), v, κ, p, ξ₁)
```
then
```
u_T(ξ₁, (p, κ, γ₀)) = γ₀ * P(ξ₁, (p, κ)) + E(ξ₁, (p, κ)) * g(ξ₁)
```
where
```
abs(g(ξ₁)) <= C₂ * exp(-real(c) * ξ₁^2) * ξ₁^(-2 / σ + 2σ * v + v + d - 2)
```
and `real(c) = κ * ϵ / (1 + ϵ)^2`.

## Enclosing `u_T_dγ(ξ₁, (p, κ, γ_ball))`
If we let
```
C₃ = C_fix_point_dγ(abs(γ₀), v, κ, p, ξ₁)
```
then
```
u_T_dγ(ξ₁, (p, κ, γ_ball)) = P(ξ₁, (p, κ)) + h(ξ₁)
```
where
```
abs(h(ξ₁)) <= C₃ * XXX
```
- **TODO**: Finish this

"""
function check_existence_fixed_point(
    κ::Arb,
    μ::Arb,
    p::AbstractGLParams{Arb},
    ξ₁::Arb,
    v::Arb;
    non_rigorous::Bool = false,
)
    # Compute u_0(ξ₁, (p, κ, μ))
    u_0 = let
        sol = if non_rigorous
            solution_zero_float(μ, κ, ξ₁, p)
        else
            solution_zero_capd(μ, κ, ξ₁, p)
        end
        Acb(sol[1], sol[2])
    end

    γ = solution_infinity_γ(κ, μ, p, ξ₁, v, u_0)

    return isfinite(γ), γ
end
