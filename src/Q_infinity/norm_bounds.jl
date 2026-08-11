"""
    NormBounds(γ, κ, ϵ, ξ₁, v, Λ, C::FunctionBounds; include_dκ = false, include_dϵ = false)

Contains bounds for norms of `Q` and its derivatives.

It contains bounds for the norms of

- `Q`
- `Q_dξ`
- `Q_dξ_dξ`
- `Q_dξ_dξ_dξ`

If either `include_dκ = true` or `include_dϵ = true` it also include bounds for:

- `Q_dγ_real`
- `Q_dγ_real_dξ`
- `Q_dγ_imag`
- `Q_dγ_imag_dξ`

If `include_dκ = true` it also include bounds for:

- `Q_dκ`
- `Q_dκ_dξ`

If `include_dϵ = true` it also include bounds for:

- `Q_dϵ`
- `Q_dϵ_dξ`

It checks all the conditions on the parameters that these lemmas
assume. If the conditions for some of the lemmas are not satisfied,
then it will return indeterminate values for those bounds. When using
this struct the bounds can therefore safely be assume to hold.

The bounds are based on the following lemmas:

- `Q`: Lemma REF(lemma:fixed-point-bounds) and Proposition
  REF(prop:fixed-point)
- `Q_dξ`, `Q_dξ_dξ` and `Q_dξ_dξ_dξ`: Lemma REF(lemma:norm-dQ)
- `Q_dγ`: Lemma REF(lemma:norm-Q-dgamma)
- `Q_dκ`: Lemma REF(lemma:norm-Q-dkappa)
- `Q_dϵ`: Lemma REF(lemma:norm-Q-depsilon)
- `Q_dγ_dξ`: Lemma REF(lemma:norm-Q-dgamma-dxi)
- `Q_dκ_dξ`: Lemma REF(lemma:norm-Q-dkappa-dxi)
- `Q_dϵ_dξ`: Lemma REF(lemma:norm-Q-depsilon-dxi)
"""
struct NormBounds
    # Always included
    Q::Arb
    Q_dξ::Arb
    Q_dξ_dξ::Arb
    Q_dξ_dξ_dξ::Arb
    # Included when either include_dκ = true or include_dϵ = true (otherwise indeterminate)
    Q_dγ_real::Arb
    Q_dγ_real_dξ::Arb
    Q_dγ_imag::Arb
    Q_dγ_imag_dξ::Arb
    # Included when include_dκ = true (otherwise indeterminate)
    Q_dκ::Arb
    Q_dκ_dξ::Arb
    # Included when include_dϵ = true (otherwise indeterminate)
    Q_dϵ::Arb
    Q_dϵ_dξ::Arb

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
    )
end

function NormBounds(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds;
    include_dκ::Bool = false,
    include_dϵ::Bool = false,
)
    (; d, σ) = Λ

    norms = NormBounds()

    # This is an implicit requirement in the paper, the norm is only
    # defined for v >= 0.
    @assert v >= 0

    # These are the requirements for Lemma
    # REF(lemma:fixed-point-bounds)
    @assert ξ₁ > 1
    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0
    @assert (2σ + 1) * v - 2 < 0

    # Only check requirements for derivatives if they are computed.
    if include_dκ || include_dϵ
        # The remaining Lemmas mostly have exactly the same requirements
        # as REF(lemma:fixed-point-bounds).
        # For Lemmas REF(lemma:norm-Q-dkappa) and REF(lemma:norm-Q-depsilon)
        # this is a requirement
        @assert v > 0
        # For Lemma REF(lemma:norm-Q-dkappa-dxi) this is a requirement
        @assert ξ₁ > exp(Arb(1))
    end

    norms.Q[] = norm_bound_Q(γ, κ, ϵ, ξ₁, v, Λ, C, CI)
    norms.Q_dξ[] = norm_bound_Q_dξ(γ, κ, ϵ, ξ₁, v, Λ, C, CI, norms)
    norms.Q_dξ_dξ[] = norm_bound_Q_dξ_dξ(γ, κ, ϵ, ξ₁, v, Λ, C, CI, norms)
    norms.Q_dξ_dξ_dξ[] = norm_bound_Q_dξ_dξ_dξ(γ, κ, ϵ, ξ₁, v, Λ, C, CI, norms)

    if include_dκ || include_dϵ
        # Note that the bounds are the same for real and imaginary
        # parts
        norms.Q_dγ_real[] = norm_bound_Q_dγ(γ, κ, ϵ, ξ₁, v, Λ, C, CI, norms)
        norms.Q_dγ_imag[] = norms.Q_dγ_real
        norms.Q_dγ_real_dξ[] = norm_bound_Q_dγ_dξ(γ, κ, ϵ, ξ₁, v, Λ, C, CI, norms)
        norms.Q_dγ_imag_dξ[] = norms.Q_dγ_real_dξ
    end

    if include_dκ
        norms.Q_dκ[] = norm_bound_Q_dκ(γ, κ, ϵ, ξ₁, v, Λ, C, CI, norms)
        norms.Q_dκ_dξ[] = norm_bound_Q_dκ_dξ(γ, κ, ϵ, ξ₁, v, Λ, C, CI, norms)
    end

    if include_dϵ
        norms.Q_dϵ[] = norm_bound_Q_dϵ(γ, κ, ϵ, ξ₁, v, Λ, C, CI, norms)
        norms.Q_dϵ_dξ[] = norm_bound_Q_dϵ_dξ(γ, κ, ϵ, ξ₁, v, Λ, C, CI, norms)
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
    norm_bound_Q(γ, κ, ϵ, ξ₁, v, Λ, C, CI)

Compute a bound for the norm of `Q` using the fixed point theorem in
Proposition REF(prop:fixed-point). For this we need to find `ρ`
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

From Proposition REF(prop:fixed-point) we have that the norm of `Q`
is bounded by `ρ`.
"""
function norm_bound_Q(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
)
    (; d, σ) = Λ
    # The requirements on the parameters for Proposition
    # REF(prop:fixed-point) are the same as for Lemma
    # REF(lemma:fixed-point-bounds). These are in turn the same as for
    # Lemma REF(lemma:I_P-I_E) and are hence checked when constructing
    # CI::IBounds.

    # In this case the solution to the ODE is exactly zero.
    iszero(γ) && return Arb(0)

    C_T = CGL.C_T(κ, ϵ, ξ₁, v, Λ, C, CI)
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
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = Λ
    return C.P_dξ * abs(γ) * ξ₁^(-v - 1) +
           C_Q_dξ(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q^(2σ + 1) * ξ₁^(2σ * v - 1)
end

function norm_bound_Q_dξ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = Λ
    return C.P_dξ_dξ * abs(γ) * ξ₁^(-v - 2) +
           (
               C_Q_dξ_dξ_1(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * ξ₁^(-1) +
               C_Q_dξ_dξ_2(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q_dξ
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
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = Λ
    return C.P_dξ_dξ_dξ * abs(γ) * ξ₁^(-v - 3) +
           (
               C_Q_dξ_dξ_dξ_1(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q^2 +
               C_Q_dξ_dξ_dξ_2(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dξ * ξ₁^(-1) +
               C_Q_dξ_dξ_dξ_3(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q_dξ^2 +
               C_Q_dξ_dξ_dξ_4(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dξ_dξ
           ) *
           norms.Q^(2σ - 1) *
           ξ₁^(2σ * v - 1)
end

# This handles both dγ_real and dγ_imag
function norm_bound_Q_dγ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = Λ
    num = C.P * ξ₁^-v
    den = (1 - (2σ + 1) * C_T(κ, ϵ, ξ₁, v, Λ, C, CI) * ξ₁^(-2 + 2σ * v) * norms.Q^2σ)

    return Arblib.ispositive(den) ? num / den : indeterminate(num)
end

# This handles both dγ_real and dγ_imag. It uses norms.Q_dγ_real, but
# this is always the same as norms.Q_dγ_imag.
function norm_bound_Q_dγ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = Λ
    return C.P_dξ * ξ₁^(-v - 1) +
           (2σ + 1) *
           C_Q_dξ(κ, ϵ, ξ₁, v, Λ, C, CI) *
           norms.Q^2σ *
           norms.Q_dγ_real *
           ξ₁^(2Λ.σ * v - 1)
end

function norm_bound_Q_dκ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = Λ
    num = (
        C_Q_dκ_1(κ, ϵ, ξ₁, v, Λ, C, CI) * abs(γ) +
        (
            C_Q_dκ_2(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q^2 +
            C_Q_dκ_3(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dξ +
            C_Q_dκ_4(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q_dξ^2 +
            C_Q_dκ_5(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dξ_dξ
        ) * norms.Q^(2σ - 1)
    )
    den = (1 - C_Q_dκ_6(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q^2σ)

    Arblib.ispositive(den) ? num / den : indeterminate(num)
end

function norm_bound_Q_dκ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = Λ
    return C.P_dκ * abs(γ) * log(ξ₁) * ξ₁^(-v - 1) +
           (
        C_Q_dξ_dκ_1(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q^2 +
        C_Q_dξ_dκ_2(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dκ +
        C_Q_dξ_dκ_3(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dξ +
        C_Q_dξ_dκ_4(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q_dξ^2 +
        C_Q_dξ_dκ_5(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dξ_dξ
    ) * norms.Q^(2σ - 1)
end

function norm_bound_Q_dϵ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = Λ
    num = (
        C_Q_dϵ_1(κ, ϵ, ξ₁, v, Λ, C, CI) * abs(γ) +
        (
            C_Q_dϵ_2(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q^2 +
            C_Q_dϵ_3(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dξ +
            C_Q_dϵ_4(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q_dξ^2 +
            C_Q_dϵ_5(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dξ_dξ
        ) * norms.Q^(2σ - 1)
    )
    den = (1 - C_Q_dϵ_6(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q^2σ)

    Arblib.ispositive(den) ? num / den : indeterminate(num)
end

function norm_bound_Q_dϵ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds,
    CI::IBounds,
    norms::NormBounds,
)
    (; σ) = Λ
    return C.P_dϵ * abs(γ) * ξ₁^(-v - 1) +
           (
        C_Q_dξ_dϵ_1(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q^2 +
        C_Q_dξ_dϵ_2(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dϵ +
        C_Q_dξ_dϵ_3(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dξ +
        C_Q_dξ_dϵ_4(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q_dξ^2 +
        C_Q_dξ_dϵ_5(κ, ϵ, ξ₁, v, Λ, C, CI) * norms.Q * norms.Q_dξ_dξ
    ) * norms.Q^(2σ - 1)
end
