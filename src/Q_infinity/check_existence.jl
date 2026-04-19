function C_T1(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}, C::FunctionBounds)
    (; d, σ) = λ
    @assert (2σ + 1) * v < 2 + 2 / σ - d
    @assert 2 / d < σ

    return C.P * C.J_E / abs((2σ + 1) * v - 2) +
           C.E * C.J_P / abs((2σ + 1) * v - 2 / σ + d - 2)
end

function M(σ::Arb)
    if isone(σ)
        return sqrt(Arb(2)) / (4 - 2sqrt(Arb(2))) + 1
    elseif isequal(σ, Arb("2.3"))
        # Precomputed value
        return Arb("[3.38307764695680109166 +/- 7.83e-21]")
    else
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

        return max(M1, M2, M3, M4)
    end
end

"""
    Q_infinity_fixed_point(γ, κ, ϵ, ξ₁, v, λ, C)

To apply the fixed point theorem in Proposition
REF(prop:Q-fixed-point) we need to find `ρ` satisfying the
inequality

```
C_P * abs(γ) * ξ₁^-v + C_T1 * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) <= ρ
```

and

```
2C_T2 * ρ^2σ * ξ₁^(-2 + 2σ * v) < 1
```

The second inequality gives us a direct upper bound for `ρ`, we take
`ρ_bound` to be a value slightly less than this upper bound.

For the first inequality we note that we can take an upper bound of
`abs(γ)`, if the inequality is satisfied for this upper bound then it
is automatically satisfied for `abs(γ)`. Let `r_1` be an upper bound
of `abs(γ)`. If `r_1 = 0` then `ρ = 0` satisfies the inequality. If
`r_1 > 0`, we consider the function

```
f(ρ) = C_P * abs(γ) * ξ₁^-v + C_T1 * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) - ρ
```

We will show that this has a unique root on the interval ``0 = ρ <
ρ_bound``. The expression is always positive at `ρ = 0`, so to the
right of the root `f(ρ)` will be negative and `ρ` hence satisfies the
inequality. The zero itself is the smallest possible `ρ` satisfying
the inequality.

To prove that there is a unique root on the interval ``0 < ρ <
ρ_bound`` we note that

```
f'(ρ) = (2σ + 1) * C_T1 * ξ₁^(-2 + 2σ * v) * ρ^2σ - 1
```

has a unique root for `ρ > 0`. It follows that that `f(ρ)` has a
unique critical point. If `f(ρ_bound)` is negative it then follows
that `f(ρ)` has a unique root on the interval.

From Proposition REF(prop:Q-fixed-point) we have that the norm of `Q`
is bounded by `ρ`.
"""
function Q_infinity_fixed_point(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; σ) = λ
    @assert v > 0 # Required for the below bounds to be valid

    # In this case the solution to the ODE is exactly zero.
    iszero(γ) && return Arb(0)

    C_T_1 = C_T1(κ, ϵ, ξ₁, v, λ, C)
    C_T_2 = M(σ) * C_T_1

    # Upper bound for ρ from second inequality. We take a value that
    # is strictly lower than this (by eps(Arb)), so that we know that
    # the strict inequality is satisfied.
    ρ_bound = lbound((2C_T_2 * ξ₁^(-2 + 2σ * v))^(-1 / 2σ) - eps(Arb))

    isfinite(ρ_bound) || return indeterminate(Arb)

    # Precompute constants
    r_1 = Arblib.abs_ubound(Arb, γ)
    w_1 = C.P * r_1 * ξ₁^-v
    w_2 = C_T_1 * ξ₁^(-2 + 2σ * v)
    f(ρ) = w_1 + w_2 * ρ^(2σ + 1) - ρ

    # Since γ is non-zero at this point we should always have w_1 > 0
    # and hence f(0) should be positive.
    @assert Arblib.ispositive(f(Arb(0)))
    # Check that f is negative at the right endpoint.
    Arblib.isnegative(f(Arb(ρ_bound))) || return indeterminate(Arb)

    # We now know that there is a unique root on the interval 0 < ρ <
    # ρ_bound.

    # We first get a rough enclosure using bisection and then refine
    # it using interval Newton.
    ρ_initial = ArbExtras.refine_root_bisection(f, Arf(0), ρ_bound, rtol = Arb(1e-3))

    return ArbExtras.refine_root(f, Arb(ρ_initial), strict = false)
end
