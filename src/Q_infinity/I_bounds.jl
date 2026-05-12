"""
    IBounds(κ, ϵ, ξ₁, Λ, C::FunctionBounds; include_dκ = false, include_dϵ = false)

Contains the constants involved in asymptotic bounds for the functions
`I_P` and `I_E` as well as their derivatives.

For `I_P`, `I_E` as well as their derivatives w.r.t. `ξ` it contains
the following constants:

- `C_I_P`
- `C_I_P_1_1`
- `C_I_P_1_2`
- `C_I_P_2_1`
- `C_I_P_2_2`
- `C_I_P_2_3`
- `C_I_P_2_4`
- `C_I_P_dξ`
- `C_I_E`
- `C_I_E_dξ`

If `include_dκ = true` it also includes the following constants
related to the derivatives of `I_P` and `I_E` w.r.t. `κ`:

- `C_I_P_dκ_1_1`
- `C_I_P_dκ_1_2`
- `C_I_P_dκ_1_3`
- `C_I_P_dκ_1_4`
- `C_I_E_dκ`

If `include_dκ = true` it also includes the following constants
related to the derivatives of `I_P` and `I_E` w.r.t. `ϵ`:

- `C_I_P_dϵ_1_1`
- `C_I_P_dϵ_1_2`
- `C_I_P_dϵ_1_3`
- `C_I_P_dϵ_1_4`
- `C_I_E_dϵ`

It checks all the conditions required for the bounds involving these
constants to be valid. If any of the checks fails it throws an error.
When using this struct the bounds can therefore safely be assume to
hold.

The bounds are based on the following lemmas:

- `C_I_P` and `C_I_E`: Lemma REF(lemma:I_P-I_E)
- `C_I_P_dξ`, `C_I_E_dξ`: Lemma REF(lemma:I_P-I_E-derivatives)
- `C_I_P_1_n` and `C_I_P_2_n`: Lemma REF(lemma:I-P-refined)
- `C_I_P_dκ_1_n` and `C_I_P_dϵ_1_n`: Lemma REF(lemma:I-P-dkappa-depsilon-1)
- `C_I_E_dκ` and `C_I_E_dϵ`: Lemma REF(lemma:I-E-bounds)
"""
struct IBounds
    # Always included
    I_P::Arb
    I_P_1_1::Arb
    I_P_1_2::Arb
    I_P_2_1::Arb
    I_P_2_2::Arb
    I_P_2_3::Arb
    I_P_2_4::Arb
    I_P_dξ::Arb
    I_E::Arb
    I_E_dξ::Arb
    # Included when include_dκ = true (otherwise indeterminate)
    I_P_dκ_1_1::Arb
    I_P_dκ_1_2::Arb
    I_P_dκ_1_3::Arb
    I_P_dκ_1_4::Arb
    I_E_dκ::Arb
    # Included when include_dϵ = true (otherwise indeterminate)
    I_P_dϵ_1_1::Arb
    I_P_dϵ_1_2::Arb
    I_P_dϵ_1_3::Arb
    I_P_dϵ_1_4::Arb
    I_E_dϵ::Arb

    IBounds() = new(
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

function IBounds(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds;
    include_dκ::Bool = false,
    include_dϵ::Bool = false,
)
    (; σ, d) = Λ

    CI = IBounds()

    # These cover all requirements of Lemmas REF(lemma:I_P-I_E),
    # REF(lemma:I_P-I_E-derivatives), REF(lemma:I-P-refined),
    # REF(lemma:I-E-bounds), REF(lemma:I-P-dkappa-1) and
    # REF(lemma:I-P-depsilon-1).
    @assert ξ₁ > 1
    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0
    @assert (2σ + 1) * v - 2 < 0

    CI.I_E[] = C_I_E(κ, ϵ, ξ₁, v, Λ, C)

    CI.I_P[] = C_I_P(κ, ϵ, ξ₁, v, Λ, C)
    CI.I_P_1_1[] = C_I_P_1_1(κ, ϵ, ξ₁, v, Λ, C)
    CI.I_P_1_2[] = C_I_P_1_2(κ, ϵ, ξ₁, v, Λ, C)
    CI.I_P_2_1[] = C_I_P_2_1(κ, ϵ, ξ₁, v, Λ, C)
    CI.I_P_2_2[] = C_I_P_2_2(κ, ϵ, ξ₁, v, Λ, C)
    CI.I_P_2_3[] = C_I_P_2_3(κ, ϵ, ξ₁, v, Λ, C)
    CI.I_P_2_4[] = C_I_P_2_4(κ, ϵ, ξ₁, v, Λ, C)

    CI.I_E_dξ[] = C_I_E_dξ(κ, ϵ, ξ₁, v, Λ, C)
    CI.I_P_dξ[] = C_I_P_dξ(κ, ϵ, ξ₁, v, Λ, C)

    if include_dκ
        CI.I_E_dκ[] = C_I_E_dκ(κ, ϵ, ξ₁, v, Λ, C)

        CI.I_P_dκ_1_1[] = C_I_P_dκ_1_1(κ, ϵ, ξ₁, v, Λ, C)
        CI.I_P_dκ_1_2[] = C_I_P_dκ_1_2(κ, ϵ, ξ₁, v, Λ, C)
        CI.I_P_dκ_1_3[] = C_I_P_dκ_1_3(κ, ϵ, ξ₁, v, Λ, C)
        CI.I_P_dκ_1_4[] = C_I_P_dκ_1_4(κ, ϵ, ξ₁, v, Λ, C)
    end

    if include_dϵ
        CI.I_E_dϵ[] = C_I_E_dϵ(κ, ϵ, ξ₁, v, Λ, C)

        CI.I_P_dϵ_1_1[] = C_I_P_dϵ_1_1(κ, ϵ, ξ₁, v, Λ, C)
        CI.I_P_dϵ_1_2[] = C_I_P_dϵ_1_2(κ, ϵ, ξ₁, v, Λ, C)
        CI.I_P_dϵ_1_3[] = C_I_P_dϵ_1_3(κ, ϵ, ξ₁, v, Λ, C)
        CI.I_P_dϵ_1_4[] = C_I_P_dϵ_1_4(κ, ϵ, ξ₁, v, Λ, C)
    end

    return CI
end

function C_I_E(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ) = Λ
    return C.J_E / abs((2σ + 1) * v - 2)
end

function C_I_P(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    bound = C.J_P / abs((2σ + 1) * v - 2 / σ + d - 2)
    if real(c) > 0
        bound2 = C.J_P / 2real(c) * ξ₁^-2
        # We return the best of the two bounds.
        return min(bound, bound2)
    else
        return bound
    end
end

function C_I_P_1_1(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return abs(B_W(κ, ϵ, Λ) / 2c) * (
        C.P +
        C.P_dξ / abs((2Λ.σ + 1) * v - 2 / Λ.σ + Λ.d - 4) +
        abs(d - 2) * C.P / abs((2Λ.σ + 1) * v - 2 / Λ.σ + Λ.d - 4)
    )
end

function C_I_P_1_2(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return abs(B_W(κ, ϵ, Λ) / 2c) * (2σ + 1) * C.P / abs((2σ + 1) * v - 2 / σ + d - 3)
end

function C_I_P_2_1(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return abs(B_W(κ, ϵ, Λ) / 2c) * C.P +
           abs(B_W(κ, ϵ, Λ) / 4c^2) *
           (
               C.P_dξ +
               abs(d - 2) * C.P +
               (C.P_dξ_dξ + abs(2d - 5) * C.P_dξ + abs((d - 2) * (d - 4)) * C.P) /
               abs((2Λ.σ + 1) * v - 2 / Λ.σ + Λ.d - 6)
           ) *
           ξ₁^-2
end

function C_I_P_2_2(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return abs(B_W(κ, ϵ, Λ) / 4c^2) *
           (2σ + 1) *
           (C.P + (2C.P_dξ + abs(2d - 5) * C.P) / abs((2Λ.σ + 1) * v - 2 / Λ.σ + Λ.d - 5))
end

function C_I_P_2_3(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return abs(B_W(κ, ϵ, Λ) / 4c^2) * (2σ + 1) * 2σ * C.P /
           abs((2Λ.σ + 1) * v - 2 / Λ.σ + Λ.d - 4)
end

function C_I_P_2_4(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return abs(B_W(κ, ϵ, Λ) / 4c^2) * (2σ + 1) * C.P /
           abs((2Λ.σ + 1) * v - 2 / Λ.σ + Λ.d - 4)
end

function C_I_E_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    return C.J_E
end

function C_I_P_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    return C.J_P
end

function C_I_E_dκ(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    return C.J_E_dκ * (inv(((2Λ.σ + 1) * v - 2)^2) + log(ξ₁) / abs((2Λ.σ + 1) * v - 2))
end

function C_I_P_dκ_1_1(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return C.D / abs(2c) +
           (
        C.D_dξ +
        d * C.D +
        (C.D_dξ_dξ + abs(2d - 1) * C.D_dξ + abs(d * (d - 2)) * C.D) /
        abs((2σ + 1) * v - 2 / σ + d - 4)
    ) / abs(2c)^2 * ξ₁^-2
end

function C_I_P_dκ_1_2(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return (2σ + 1) *
           (C.D + (2C.D_dξ + abs(2d - 1) * C.D) / abs((2σ + 1) * v - 2 / σ + d - 3)) /
           abs(2c)^2
end

function C_I_P_dκ_1_3(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return (2σ + 1) * 2σ * C.D / abs((2σ + 1) * v - 2 / σ + d - 2) / abs(2c)^2
end

function C_I_P_dκ_1_4(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return (2σ + 1) * C.D / abs((2σ + 1) * v - 2 / σ + d - 2) / abs(2c)^2
end

function C_I_E_dϵ(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    return C.J_E_dϵ / abs((2Λ.σ + 1) * v - 2)
end

function C_I_P_dϵ_1_1(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return C.H / abs(2c) +
           (
        C.H_dξ +
        d * C.H +
        (C.H_dξ_dξ + abs(2d - 1) * C.H_dξ + abs(d * (d - 2)) * C.H) /
        abs((2σ + 1) * v - 2 / σ + d - 4)
    ) / abs(2c)^2 * ξ₁^-2
end

function C_I_P_dϵ_1_2(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return (2σ + 1) *
           (C.H + (2C.H_dξ + abs(2d - 1) * C.H) / abs((2σ + 1) * v - 2 / σ + d - 3)) /
           abs(2c)^2
end

function C_I_P_dϵ_1_3(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return (2σ + 1) * 2σ * C.H / abs((2σ + 1) * v - 2 / σ + d - 2) / abs(2c)^2
end

function C_I_P_dϵ_1_4(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, Λ::CGLParams{Arb}, C::FunctionBounds)
    (; σ, d) = Λ
    c = _c(κ, ϵ, Λ)
    return (2σ + 1) * C.H / abs((2σ + 1) * v - 2 / σ + d - 2) / abs(2c)^2
end
