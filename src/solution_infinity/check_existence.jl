function C_T1(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; d, σ) = λ

    @assert (2σ + 1) * v < 2 + 2 / σ - d
    @assert 2 / d < σ

    return C_P(κ, λ, ξ₁) * C_J_E(κ, ξ₁, λ) / abs((2σ + 1) * v - 2) +
           C_E(κ, λ, ξ₁) * C_J_P(κ, ξ₁, λ) / abs((2σ + 1) * v - 2 / σ + d - 2)
end

function C_T2(κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    M = if isone(λ.σ)
        sqrt(Arb(2)) / (4 - 2sqrt(Arb(2))) + 1
    else
        (; σ) = λ

        t₀ = Arb(1 - 1e-3)
        t₁ = Arb(1 + 1e-3)

        g(t) = (1 - abspow(t, 2σ)) / ((1 - t) * (1 + abspow(t, 2σ))) + 1

        # Bound for 0 <= t <= t₀
        M1 = ArbExtras.maximum_enclosure(g, Arf(0), ubound(t₀))
        # Bound for t₀ <= t <= t₁
        M2 = (1 - t₁^2σ) / ((1 - t₁) * (1 + t₀^2σ)) + 1
        # Bound for t₁ <= t <= 2
        M3 = ArbExtras.maximum_enclosure(g, lbound(t₁), Arf(2))
        # Bound for t >= 2
        M4 = Arb(2)

        max(M1, M2, M3, M4)
    end

    return M * C_T1(κ, ξ₁, v, λ)
end

"""
    solution_infinity_fixed_point(γ, κ, ξ₁, v, λ::CGLParams)

Consider the fixed point problem given by
[`fpp_infinity_complex`](@ref). This function computes `ρ_l, ρ_u` such
that there exists a unique fixed point in the ball of radius `ρ_u` and
it is contained in the ball of radius `ρ_l`.

We have a unique fixed points in a ball of radius `ρ` if
```
C_P * r1 * ξ₁^-v + C_T1 * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) <= ρ
```
and
```
2C_T2 * ρ^2σ * ξ₁^(-2 + 2σ * v) < 1
```
where `r1 = abs(γ)`.

The second inequality gives us a direct upper bound for `ρ`. For the
first inequality we find the zeros of
```
C_P * r1 * ξ₁^-v + C_T1 * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) - ρ
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

    @assert v > 0 # Required for the below bounds to be valid

    r1 = abs(γ)

    CP = C_P(κ, λ, ξ₁)
    CT1 = C_T1(κ, ξ₁, v, λ)
    CT2 = C_T2(κ, ξ₁, v, λ)

    # Upper from second inequality
    ρ_bound = (2CT2 * ξ₁^(-2 + 2σ * v))^(-1 / 2σ)

    f(ρ) = CP * r1 * ξ₁^-v + CT1 * ξ₁^(-2 + 2σ * v) * abspow(ρ, 2σ + 1) - ρ

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
