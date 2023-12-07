"""
    solution_infinity_fixed_point(γ, κ, ξ₁, v, λ::CGLParams)

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
    λ::CGLParams{Arb};
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
