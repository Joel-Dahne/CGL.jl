"""
    NormBounds(γ, κ, ϵ, ξ₁, v, λ, C::FunctionBounds; include_dκ = false, include_dϵ = false)

Contains bounds for norms of `Q` and its derivatives. It also contains
specific bounds for `abs(Q)^2σ * Q` and its derivatives.

It contains bounds for the norms of

- `Q`
- `Q_dξ`
- `Q_dξ_dξ`
- `Q_dξ_dξ_dξ`
- `abs(Q)^2σ * Q`
- `abs(Q)^2σ * Q` differentiated one time w.r.t. `ξ`
- `abs(Q)^2σ * Q` differentiated two times w.r.t. `ξ`
- `abs(Q)^2σ * Q` differentiated three times w.r.t. `ξ`

If either `include_dκ = true` or `include_dϵ = true` it also include bounds for:

- `Q_dγ`
- `Q_dγ_dξ`
- `abs(Q)^2σ * Q` differentiated w.r.t. `γ`
- `abs(Q)^2σ * Q` differentiated w.r.t. `γ` and `ξ`

If `include_dκ = true` it also include bounds for:

- `Q_dκ`
- `Q_dκ_dξ`
- `abs(Q)^2σ * Q` differentiated w.r.t. `κ`
- `abs(Q)^2σ * Q` differentiated w.r.t. `κ` and `ξ`

If `include_dϵ = true` it also include bounds for:

- `Q_dϵ`
- `Q_dϵ_dξ`
- `abs(Q)^2σ * Q` differentiated w.r.t. `ϵ`
- `abs(Q)^2σ * Q` differentiated w.r.t. `ϵ` and `ξ`

It checks all the conditions on the parameters that these lemmas
assume. If the conditions for some of the lemmas are not satisfied,
then it will return indeterminate values for those bounds. When using
this struct the bounds can therefore safely be assume to hold.

TODO: Add documentation about which lemmas are used.
"""
struct NormBounds
    # Always included
    Q::Arb
    Q_dξ::Arb
    Q_dξ_dξ::Arb
    Q_dξ_dξ_dξ::Arb
    Q2σQ::Arb
    Q2σQ_dξ::Arb
    Q2σQ_dξ_dξ::Arb
    Q2σQ_dξ_dξ_dξ::Arb
    # Included when either include_dκ = true or include_dϵ = true (otherwise indeterminate)
    Q_dγ::Arb
    Q_dγ_dξ::Arb
    Q2σQ_dγ::Arb
    Q2σQ_dγ_dξ::Arb
    # Included when include_dκ = true (otherwise indeterminate)
    Q_dκ::Arb
    Q_dκ_dξ::Arb
    Q2σQ_dκ::Arb
    Q2σQ_dκ_dξ::Arb
    # Included when include_dϵ = true (otherwise indeterminate)
    Q_dϵ::Arb
    Q_dϵ_dξ::Arb
    Q2σQ_dϵ::Arb
    Q2σQ_dϵ_dξ::Arb

    NormBounds() = new(
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
    )
end

function NormBounds(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds;
    include_dκ::Bool = false,
    include_dϵ::Bool = false,
)
    norms = NormBounds()

    (; σ) = λ

    norms.Q[] = norm_bound_Q(γ, κ, ϵ, ξ₁, v, λ, C, CI)
    norms.Q_dξ[] = norm_bound_Q_dξ(γ, κ, ϵ, ξ₁, v, λ, C, CI, norms)
    norms.Q_dξ_dξ[] = norm_bound_Q_dξ_dξ(γ, κ, ϵ, ξ₁, v, λ, C, CI, norms)
    norms.Q_dξ_dξ_dξ[] = norm_bound_Q_dξ_dξ_dξ(γ, κ, ϵ, ξ₁, v, λ, C, CI, norms)

    # Norms of abs(Q)^2σ * Q and its derivatives
    # TODO: Bounds this as abs(Q)^(2σ + 1)?
    norm_Q_series =
        ArbSeries((norms.Q, norms.Q_dξ, norms.Q_dξ_dξ / 2, norms.Q_dξ_dξ_dξ / 6))
    norm_Q2σQ_series = norm_Q_series^2σ * norm_Q_series

    norms.Q2σQ[] = norm_Q2σQ_series[0]
    norms.Q2σQ_dξ[] = norm_Q2σQ_series[1]
    norms.Q2σQ_dξ_dξ[] = 2norm_Q2σQ_series[2]
    norms.Q2σQ_dξ_dξ_dξ[] = 6norm_Q2σQ_series[3]

    if include_dκ || include_dϵ
        norms.Q_dγ[] = norm_bound_Q_dγ(γ, κ, ϵ, ξ₁, v, λ, C, CI, norms)
        norms.Q_dγ_dξ[] = norm_bound_Q_dγ_dξ(γ, κ, ϵ, ξ₁, v, λ, C, CI, norms)

        # Norms of abs(Q)^2σ * Q differentiated w.r.t. γ
        norms.Q2σQ_dγ[] = (2σ + 1) * norms.Q^2σ * norms.Q_dγ
        norms.Q2σQ_dγ_dξ[] =
            (2σ + 1) *
            norms.Q^(2σ - 1) *
            (2σ * norms.Q_dξ * norms.Q_dγ + norms.Q * norms.Q_dγ_dξ)
    end

    if include_dκ
        norms.Q_dκ[] = norm_bound_Q_dκ(γ, κ, ϵ, ξ₁, v, λ, C, CI, norms)
        norms.Q_dκ_dξ[] = norm_bound_Q_dκ_dξ(γ, κ, ϵ, ξ₁, v, λ, C, CI, norms)

        # Norms of abs(Q)^2σ * Q differentiated w.r.t. κ
        norms.Q2σQ_dκ[] = (2σ + 1) * norms.Q^2σ * norms.Q_dκ
        norms.Q2σQ_dκ_dξ[] =
            (2σ + 1) *
            norms.Q^(2σ - 1) *
            (2σ * norms.Q_dξ * norms.Q_dκ + norms.Q * norms.Q_dκ_dξ)
    end

    if include_dϵ
        norms.Q_dϵ[] = norm_bound_Q_dϵ(γ, κ, ϵ, ξ₁, v, λ, C, CI, norms)
        norms.Q_dϵ_dξ[] = norm_bound_Q_dϵ_dξ(γ, κ, ϵ, ξ₁, v, λ, C, CI, norms)

        # Norms of abs(Q)^2σ * Q differentiated w.r.t. ϵ
        norms.Q2σQ_dϵ[] = (2σ + 1) * norms.Q^2σ * norms.Q_dϵ
        norms.Q2σQ_dϵ_dξ[] =
            (2σ + 1) *
            norms.Q^(2σ - 1) *
            (2σ * norms.Q_dξ * norms.Q_dϵ + norms.Q * norms.Q_dϵ_dξ)
    end

    return norms
end

"""
    M(σ::Arb)

Compute an enclosure of `M_σ` from Lemma REF(lemma:M).
"""
function M(σ::Arb)
    if isone(σ)
        return sqrt(Arb(2)) / (4 - 2sqrt(Arb(2))) + 1
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
    norm_bound_Q(γ, κ, ϵ, ξ₁, v, λ, C, CI)

Compute a bound for the norm of `Q` using the fixed point theorem in
Proposition REF(prop:Q-fixed-point). For this we need to find `ρ`
satisfying the inequality

```
C_P * abs(γ) * ξ₁^-v + C_T * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) <= ρ
```

and

```
2M(σ) * C_T * ρ^2σ * ξ₁^(-2 + 2σ * v) < 1
```

The second inequality gives us a direct upper bound for `ρ`, we take
`ρ_bound` to be a value slightly less than this upper bound.

For the first inequality we note that we can take an upper bound of
`abs(γ)`, if the inequality is satisfied for this upper bound then it
is automatically satisfied for `abs(γ)`. Let `r_1` be an upper bound
of `abs(γ)`. If `r_1 = 0` then `ρ = 0` satisfies the inequality. If
`r_1 > 0`, we consider the function

```
f(ρ) = C_P * abs(γ) * ξ₁^-v + C_T * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) - ρ
```

We will show that this has a unique root on the interval ``0 = ρ <
ρ_bound``. The expression is always positive at `ρ = 0`, so to the
right of the root `f(ρ)` will be negative and `ρ` hence satisfies the
inequality. The zero itself is the smallest possible `ρ` satisfying
the inequality.

To prove that there is a unique root on the interval ``0 < ρ <
ρ_bound`` we note that

```
f'(ρ) = (2σ + 1) * C_T * ξ₁^(-2 + 2σ * v) * ρ^2σ - 1
```

has a unique root for `ρ > 0`. It follows that that `f(ρ)` has a
unique critical point. If `f(ρ_bound)` is negative it then follows
that `f(ρ)` has a unique root on the interval.

From Proposition REF(prop:Q-fixed-point) we have that the norm of `Q`
is bounded by `ρ`.
"""
function norm_bound_Q(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; d, σ) = λ
    # Check conditions for lemma
    @assert v > 0
    @assert (2σ + 1) * v < 2 + 2 / σ - d
    @assert 2 / d < σ

    # In this case the solution to the ODE is exactly zero.
    iszero(γ) && return Arb(0)

    C_T = CGL.C_T(κ, ϵ, ξ₁, v, λ, C, CI)
    M_σ = M(σ)

    # Upper bound for ρ from second inequality. We take a value that
    # is strictly lower than this (by eps(Arb)), so that we know that
    # the strict inequality is satisfied.
    ρ_bound = lbound((2M_σ * C_T * ξ₁^(-2 + 2σ * v))^(-1 / 2σ) - eps(Arb))

    # Verify that ρ_bound is positive
    ρ_bound > 0 || return indeterminate(Arb)

    # Precompute constants
    r_1 = Arblib.abs_ubound(Arb, γ)
    w_1 = C.P * r_1 * ξ₁^-v
    w_2 = C_T * ξ₁^(-2 + 2σ * v)
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

function norm_bound_Q_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = λ
    return C.P_dξ * abs(γ) * ξ₁^(-v - 1) +
           C_Q_dξ(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q^(2σ + 1) * ξ₁^(2σ * v - 1)
end

function norm_bound_Q_dξ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = λ
    return C.P_dξ_dξ * abs(γ) * ξ₁^(-v - 2) +
           (
               C_Q_dξ_dξ_1(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * ξ₁^(-1) +
               C_Q_dξ_dξ_2(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q_dξ
           ) *
           norms.Q^2σ *
           ξ₁^(2σ * v - 1)
end

function norm_bound_Q_dξ_dξ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = λ
    return C.P_dξ_dξ_dξ * abs(γ) * ξ₁^(-v - 3) +
           (
               C_Q_dξ_dξ_dξ_1(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q^2 +
               C_Q_dξ_dξ_dξ_2(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dξ * ξ₁^(-1) +
               C_Q_dξ_dξ_dξ_3(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q_dξ^2 +
               C_Q_dξ_dξ_dξ_4(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dξ_dξ
           ) *
           norms.Q^(2σ - 1) *
           ξ₁^(2σ * v - 1)
end

function norm_bound_Q_dγ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = λ
    num = C.P * ξ₁^-v
    den = (1 - (2σ + 1) * C_T(κ, ϵ, ξ₁, v, λ, C, CI) * ξ₁^(-2 + 2σ * v) * norms.Q^2σ)

    return Arblib.ispositive(den) ? num / den : indeterminate(num)
end

function norm_bound_Q_dγ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = λ
    return C.P_dξ * ξ₁^(-v - 1) +
           (2σ + 1) *
           C_Q_dξ(κ, ϵ, ξ₁, v, λ, C, CI) *
           norms.Q^2σ *
           norms.Q_dγ *
           ξ₁^(2λ.σ * v - 1)
end

function norm_bound_Q_dκ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = λ
    num = (
        C_Q_dκ_1(κ, ϵ, ξ₁, v, λ, C, CI) * abs(γ) +
        (
            C_Q_dκ_2(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q^2 +
            C_Q_dκ_3(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dξ +
            C_Q_dκ_4(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q_dξ^2 +
            C_Q_dκ_5(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dξ_dξ
        ) * norms.Q^(2σ - 1)
    )
    den = (1 - C_Q_dκ_6(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q^2σ)

    Arblib.ispositive(den) ? num / den : indeterminate(num)
end

function norm_bound_Q_dκ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = λ
    @assert ξ₁ >= ℯ
    return C.P_dκ * abs(γ) * log(ξ₁) * ξ₁^(-v - 1) +
           (
        C_Q_dξ_dκ_1(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q^2 +
        C_Q_dξ_dκ_2(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dκ +
        C_Q_dξ_dκ_3(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dξ +
        C_Q_dξ_dκ_4(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q_dξ^2 +
        C_Q_dξ_dκ_5(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dξ_dξ
    ) * norms.Q^(2σ - 1)
end

function norm_bound_Q_dϵ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = λ
    num = (
        C_Q_dϵ_1(κ, ϵ, ξ₁, v, λ, C, CI) * abs(γ) +
        (
            C_Q_dϵ_2(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q^2 +
            C_Q_dϵ_3(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dξ +
            C_Q_dϵ_4(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q_dξ^2 +
            C_Q_dϵ_5(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dξ_dξ
        ) * norms.Q^(2σ - 1)
    )
    den = (1 - C_Q_dϵ_6(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q^2σ)

    Arblib.ispositive(den) ? num / den : indeterminate(num)
end

function norm_bound_Q_dϵ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = λ
    return C.P_dϵ * abs(γ) * ξ₁^(-v - 1) +
           (
        C_Q_dξ_dϵ_1(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q^2 +
        C_Q_dξ_dϵ_2(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dϵ +
        C_Q_dξ_dϵ_3(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dξ +
        C_Q_dξ_dϵ_4(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q_dξ^2 +
        C_Q_dξ_dϵ_5(κ, ϵ, ξ₁, v, λ, C, CI) * norms.Q * norms.Q_dξ_dξ
    ) * norms.Q^(2σ - 1)
end
