"""
    check_existence_fixed_point(
        μ::Arb,
        κ::Arb,
        ξ₁::Arb,
        v::Arb,
        λ::AbstractGLParams{Arb};
        non_rigorous::Bool = false,
    )

Let `Q_0` be the solution to [`ivp_zero_complex`](@ref). This function
attempt to prove that the there exists `γ` such that the operator `T`
has a fixed point `u_T(ξ, (p, κ, γ))` satisfying that
```
Q_0(ξ₁) = u_T(ξ₁, (p, κ, γ))
```
It returns `true, γ` on success, where `γ` is a ball enclosing the
existing value of `γ`. On failure it returns false and an
indeterminate value for `γ`.

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
    μ::Arb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::AbstractGLParams{Arb};
    non_rigorous::Bool = false,
)
    # Compute Q(ξ₁)
    Q_0 = if non_rigorous
        Acb(solution_zero_float(μ, κ, ξ₁, λ)[1:2]...)
    else
        solution_zero(μ, κ, ξ₁, λ)[1]
    end

    γ = solution_infinity_γ(κ, μ, λ, ξ₁, v, Q_0)

    return isfinite(γ), γ
end

"""
    solution_infinity_fixed_point(γ, κ, ξ₁, v, λ::AbstractGLParams)

Consider the fixed point problem given by
[`fpp_infinity_complex`](@ref). This function computes `ρ_l, ρ_u` such
that there exists a unique fixed point in the ball of radius `ρ_u` and
it is contained in the ball of radius `ρ_l`.

We have a unique fixed points in a ball of radius `ρ` if
```
C_T1 * r1 * ξ₁^-v + C_T2 * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) <= ρ
```
and
```
2C_T3 * ρ^2σ * ξ₁^(-2 + 2σ * v) < 1
```
where `r1 = abs(γ)`.

The second inequality gives us a direct upper bound for `ρ`. For the
first inequality we find the zeros of
```
C_T1 * r1 * ξ₁^-v + C_T2 * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) - ρ
```
This is always positive at `ρ = 0` so it is enough to find the
smallest root to find `ρ_l`.
"""
function solution_infinity_fixed_point(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::AbstractGLParams{Arb};
    throw_on_failure::Bool = true,
)
    (; σ) = λ

    r1 = abs(γ)

    CT1, CT2 = C_T1(v, κ, λ, ξ₁)
    CT3 = C_T2(v, κ, λ, ξ₁)

    # Upper from second inequality
    ρ_bound = (2CT3 * ξ₁^(-2 + 2σ * v))^(-1 / 2σ)

    f(ρ) = CT1 * r1 * ξ₁^-v + CT2 * ξ₁^(-2 + 2σ * v) * abspow(ρ, 2σ + 1) - ρ

    # Isolate roots
    roots, flags = ArbExtras.isolate_roots(f, Arf(0), ubound(ρ_bound))

    if length(roots) == 1 && only(flags)
        ρ_l = ArbExtras.refine_root(f, Arb(only(roots)))
        ρ_u = ρ_bound
    elseif length(roots) == 2 && all(flags)
        ρ_l = ArbExtras.refine_root(f, Arb(roots[1]))
        ρ_u = min(ρ_bound, ArbExtras.refine_root(f, Arb(roots[2])))
    elseif throw_on_failure
        if isempty(roots)
            error("could not find any roots when bounding fixed point")
        elseif !all(flags)
            error("could not isolate roots when bounding fixed point")
        else
            error("found more than two roots when bounding fixed point")
        end
    else
        ρ_l = indeterminate(ρ_bound)
        ρ_u = indeterminate(ρ_bound)
    end

    return ρ_l, ρ_u
end
