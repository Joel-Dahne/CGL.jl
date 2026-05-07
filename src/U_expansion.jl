"""
    UBounds(a, b, c, ξ₁; include_da = false)

Contains the constants involved in asymptotic bounds for the confluent
hypergeometric function `U` and related derivatives.

It contains asymptotic bounds:

```
U(a, b, z₁)
U(a + 1, b + 1, z₁)
U(a + 2, b + 2, z₁)
U(a + 3, b + 3, z₁)
U(b - a, b, -z₁)
U(b - a + 1, b + 1, -z₁)
U(b - a + 2, b + 2, -z₁)
U(b - a + 3, b + 3, -z₁)
```

These bounds are based on Lemma REF(lemma:U).

It also contains asymptotic bounds for derivatives w.r.t. `z`, namely:

```
U_dz(a, b, z₁)
U_dz(a, b, z₁, 2)
U_dz(a, b, z₁, 3)
U_dz(b - a, b, -z₁)
U_dz(b - a, b, -z₁, 2)
U_dz(b - a, b, -z₁, 3)
```

These are computed based on the earlier bounds, using that the `n`th
derivative w.r.t. `z` is given by

```
U_dz(a, b, z₁, n) = (-1)^n * U(a + n, b + n, z₁) * rising(a, n)
```

If `include_da = true` it also includes bounds for:

```
U_da(a, b, z₁)
U_da(a + 1, b + 1, z₁)
U_da(a + 2, b + 2, z₁)
U_da(b - a, b, -z₁)
U_da(b - a + 1, b + 1, -z₁)
```

These bounds are based on Lemma REF(lemma:U-a).

It checks all the conditions on the parameters that these lemmas
assume. If any of these conditions are not satisfied it will return
infinite (indeterminate) bounds. When using this struct, the bounds
can therefore safely be assume to hold.
"""
struct UBounds
    # U
    U_a_b::Arb
    U_ap1_bp1::Arb
    U_ap2_bp2::Arb
    U_ap3_bp3::Arb
    U_bma_b::Arb
    U_bmap1_bp1::Arb
    U_bmap2_bp2::Arb
    U_bmap3_bp3::Arb
    # U_dz
    U_dz_a_b::Arb
    U_dz_bma_b::Arb
    # U_dz_dz
    U_dz_dz_a_b::Arb
    U_dz_dz_bma_b::Arb
    # U_dz_dz_dz
    U_dz_dz_dz_a_b::Arb
    U_dz_dz_dz_bma_b::Arb
    # U_da
    U_da_a_b::Arb
    U_da_ap1_bp1::Arb
    U_da_ap2_bp2::Arb
    U_da_bma_b::Arb
    U_da_bmap1_bp1::Arb

    UBounds() = new(
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

function UBounds(a::Acb, b::Acb, c::Acb, ξ₁::Arb; include_da::Bool = false)
    CU = UBounds()

    z₁ = c * ξ₁^2
    mz₁ = -z₁

    CU.U_a_b[] = C_U(a, b, z₁)
    CU.U_ap1_bp1[] = C_U(a + 1, b + 1, z₁)
    CU.U_ap2_bp2[] = C_U(a + 2, b + 2, z₁)
    CU.U_ap3_bp3[] = C_U(a + 3, b + 3, z₁)
    CU.U_bma_b[] = C_U(b - a, b, mz₁)
    CU.U_bmap1_bp1[] = C_U(b - a + 1, b + 1, mz₁)
    CU.U_bmap2_bp2[] = C_U(b - a + 2, b + 2, mz₁)
    CU.U_bmap3_bp3[] = C_U(b - a + 3, b + 3, mz₁)

    CU.U_dz_a_b[] = CU.U_ap1_bp1 * abs(rising(a, 1))
    CU.U_dz_dz_a_b[] = CU.U_ap2_bp2 * abs(rising(a, 2))
    CU.U_dz_dz_dz_a_b[] = CU.U_ap3_bp3 * abs(rising(a, 3))
    CU.U_dz_bma_b[] = CU.U_bmap1_bp1 * abs(rising(b - a, 1))
    CU.U_dz_dz_bma_b[] = CU.U_bmap2_bp2 * abs(rising(b - a, 2))
    CU.U_dz_dz_dz_bma_b[] = CU.U_bmap3_bp3 * abs(rising(b - a, 3))

    if include_da
        CU.U_da_a_b[] = C_U_da(a, b, z₁)
        CU.U_da_ap1_bp1[] = C_U_da(a + 1, b + 1, z₁)
        CU.U_da_ap2_bp2[] = C_U_da(a + 2, b + 2, z₁)
        CU.U_da_bma_b[] = C_U_da(b - a, b, mz₁)
        CU.U_da_bmap1_bp1[] = C_U_da(b - a + 1, b + 1, mz₁)
    end

    return CU
end

# Coefficients in asymptotic expansion

"""
    p_U!(res::Acb, k::Integer, a::Acb, b::Acb, z::Acb, tmp::Acb = Acb())

Compute [`p_U`](@ref) using mutablie arithmetic, placing the result in
`res`. The argument `tmp` is used as scratch space for the
computations.
"""
function p_U!(res::Acb, k::Integer, a::Acb, b::Acb, z::Acb, tmp::Acb = Acb())
    Arblib.sub!(tmp, a, b)
    Arblib.add!(tmp, tmp, 1)
    Arblib.rising!(res, tmp, convert(UInt, k))
    Arblib.rising!(tmp, a, convert(UInt, k))
    Arblib.mul!(res, res, tmp)
    Arblib.neg!(tmp, z)
    Arblib.pow!(tmp, tmp, k)
    Arblib.mul!(tmp, tmp, factorial(k))
    return Arblib.div!(res, res, tmp)
end

"""
    p_U(k::Integer, a, b, z)

Compute the `k`th coefficient in the asymptotic expansion of the
confluent hypergeometric function `U`. It is given by

```
rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-z)^k)
```
"""
p_U(k::Integer, a, b, z) = rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-z)^k)
p_U(k::Integer, a::Acb, b::Acb, z::Acb) = p_U!(zero(z), k, a, b, z, Acb())

"""
    p_U_p_U_da!(res::Acb, res_da::Acb, k::Integer, a::Acb, b::Acb, z::Acb, tmp1::Acb = Acb(), tmp2::Acb = Acb(), tmp3::Acb = Acb())

Compute [`p_U_p_U_da`](@ref) using mutablie arithmetic, placing the
value for `p_U` in `res` and the derivative in `res_da`. The arguments
`tmp1`, `tmp2` and `tmp3` are used as scratch space for the
computations and should not be aliased with `res` or `res_da`.
"""
function p_U_p_U_da!(
    res::Acb,
    res_da::Acb,
    k::Integer,
    a::Acb,
    b::Acb,
    z::Acb,
    tmp1::Acb = Acb(),
    tmp2::Acb = Acb(),
    tmp3::Acb = Acb(),
)
    # Compute rising(a - b + 1, k) and its derivative
    Arblib.sub!(tmp2, a, b)
    Arblib.add!(tmp2, tmp2, 1)
    Arblib.rising2!(res, tmp1, tmp2, convert(UInt, k))

    # Compute rising(a, k) and its derivative
    Arblib.rising2!(tmp2, tmp3, a, convert(UInt, k))

    # Compute the numerator for the derivative
    Arblib.mul!(res_da, res, tmp3)
    Arblib.addmul!(res_da, tmp1, tmp2)

    # Compute the numerator for the value
    Arblib.mul!(res, res, tmp2)

    # Divide by factorial(k) * (-z)^k
    Arblib.neg!(tmp1, z)
    Arblib.pow!(tmp1, tmp1, k)
    Arblib.mul!(tmp1, tmp1, factorial(k))
    return Arblib.div!(res, res, tmp1), Arblib.div!(res_da, res_da, tmp1)
end

"""
    p_U_p_U_da(k::Integer, a, b, z)

Compute both `p_U` from [`p_U`](@ref) as well as its derivative w.r.t.
`a`. These coefficient come up in the asymptotic expansion for the
derivative if `U` w.r.t. `a`. The derivative w.r.t. `a` is given by

```
(drising(a, k) * rising(a - b + 1, k) + rising(a, k) * drising(a - b + 1, k)) / (factorial(k) * (-z)^k)
```

Where we use `drising` do denote the derivative w.r.t. the first
argument.
"""
p_U_p_U_da(k::Integer, a::Acb, b::Acb, z::Acb) =
    p_U_p_U_da!(zero(z), zero(z), k, a, b, z, Acb(), Acb(), Acb())

# Bound for remainder terms in asymptotic expansion

"""
    C_R_U(n::Integer, a::Acb, b::Acb, z::Acb)

Compute the constant ``C_{R_U}`` from Lemma REF(lemma:U). See also
https://fungrim.org/entry/461a54/. It bounds the remainder term in the
asymptotic expansion of the confluent hypergeometric function `U`.
"""
function C_R_U(n::Integer, a::Acb, b::Acb, z::Acb)
    # This is a requirement of Lemma REF(lemma:U)
    abs(imag(z)) > abs(b - 2a) || abs(real(z)) > abs(b - 2a) || return indeterminate(Arb)

    σ = abs(b - 2a) / abs(z)
    ρ = abs(a^2 - a * b + b / 2) + σ * (1 + σ / 4) / (1 - σ)^2

    return abs(rising(a, n) * rising(a - b + 1, n) / factorial(n)) *
           2sqrt(1 + Arb(π) * (n // 2)) / (1 - σ) * exp(π * ρ / ((1 - σ) * abs(z)))
end

"""
    C_R_U_12(n::Integer, a::Acb, b::Acb, z::Acb)

Compute the constants ``C_{R_U}`` and ``C_{R_U}`` from Lemma
REF(lemma:U-a). The appear in the bound for the the remainder term in
the asymptotic expansion of the derivative of the confluent
hypergeometric function `U` w.r.t. `a`.
"""
function C_R_U_12(n::Integer, a::Acb, b::Acb, z::Acb)
    # These are requirements of Lemma REF(lemma:U-a)
    abs(imag(z)) > abs(b - 2a) ||
        abs(real(z)) > abs(b - 2a) ||
        return indeterminate(Arb), indeterminate(Arb)
    0 < real(a) < real(b) || return indeterminate(Arb), indeterminate(Arb)
    abs(angle(z)) < π || return indeterminate(Arb), indeterminate(Arb)
    real(a - b + n + 1) > 0 || return indeterminate(Arb), indeterminate(Arb)

    γ = angle(z)

    ρ_γ = if real(z) >= 0 # Corresponds to abs(γ) <= π / 2
        Arb(1)
    else
        # This is always >= 1, so correct even when z overlaps the
        # imaginary axis.
        inv(sin(π - abs(γ)))
    end

    C_χ =
        abs(sinpi(a - b + 1)) / π * gamma(-real(a - b)) * gamma(real(a - b + n + 1)) /
        factorial(n)

    C_χ_a =
        abs(cospi(a - b + 1)) * gamma(-real(a - b)) * gamma(real(a - b + n + 1)) /
        factorial(n) +
        abs(sinpi(a - b + 1)) / π * (
            inv(real(a - b)^2) +
            gamma(-real(a - b)) * gamma(real(a - b + n)) / factorial(n - 1)
        )

    C1 = gamma(real(a + n))

    C2 = inv(real(a + n)^2) + gamma(real(a + n + 1))

    C_R_U_1 =
        ρ_γ * C_χ / abs(gamma(a)) *
        exp(-γ * imag(a)) *
        (C1 + (abs(γ) * C1 + C2) / abs(log(abs(z))))

    C_R_U_2 = ρ_γ * C_χ_a * C1 / abs(gamma(a)) * exp(-γ * imag(a)) / abs(log(abs(z)))

    return C_R_U_1, C_R_U_2
end

"""
    C_U(a::Acb, b::Acb, z₁::Acb, n::Integer = 20)

Compute the constant `C_U` from Lemma REF(lemma:U). It satisfies that

```
abs(U(a, b, z)) <= C * abs(z^-a)
```

for `z` such that `abs(imag(z)) > abs(imag(z₁))` and `abs(z) > abs(z₁)`.

It requires that `abs(imag(z₁)) > abs(imag(b - 2a))`. If this is not
satisfied it returns an indeterminate result.
"""
function C_U(a::Acb, b::Acb, z₁::Acb, n::Integer = 20)
    isfinite(a) && isfinite(b) && isfinite(z₁) || return indeterminate(Arb)
    # This is a requirement of Lemma REF(lemma:U)
    abs(imag(z₁)) > abs(b - 2a) || abs(real(z₁)) > abs(b - 2a) || return indeterminate(Arb)

    term = zero(a)
    abs_term = zero(Arb)
    tmp = zero(a)

    S = zero(Arb)
    for k = 0:(n-1)
        Arblib.abs!(abs_term, p_U!(term, k, a, b, z₁, tmp))
        Arblib.add!(S, S, abs_term)
    end

    return S + C_R_U(n, a, b, z₁) * abs(z₁)^-n
end

"""
    C_U_da(a::Acb, b::Acb, z₁::Acb, n::Integer = 20)

Compute the constant ``C_{U,a}`` from Lemma REF(lemma:U-a). It
satisfies that

```
abs(U_da(a, b, z)) <= C_U_da * abs(log(z) * z^-a)
```

for `z` such that `abs(imag(z)) > abs(imag(z₁))` and `abs(z) > abs(z₁)`.

It requires that

```
abs(imag(z₁)) > abs(imag(b - 2a))
0 < real(a) < real(b)
abs(angle(z₁)) < π
real(a - b + n + 1) > 0
```

If any of these are not satisfied it returns an indeterminate result.
"""
function C_U_da(a::Acb, b::Acb, z₁::Acb, n::Integer = 20)
    isfinite(a) && isfinite(b) && isfinite(z₁) || return indeterminate(Arb)
    # These are requirements of Lemma REF(lemma:U-a)
    abs(imag(z₁)) > abs(b - 2a) || abs(real(z₁)) > abs(b - 2a) || return indeterminate(Arb)
    0 < real(a) < real(b) || return indeterminate(Arb)
    abs(angle(z₁)) < π || return indeterminate(Arb)
    real(a - b + n + 1) > 0 || return indeterminate(Arb)

    S1 = zero(Arb)
    S2 = zero(Arb)
    term1 = zero(a)
    term2 = zero(a)
    abs_term = zero(Arb)
    tmp1 = zero(a)
    tmp2 = zero(a)
    tmp3 = zero(a)
    for k = 0:(n-1)
        p_U_p_U_da!(term1, term2, k, a, b, z₁, tmp1, tmp2, tmp3)
        Arblib.abs!(abs_term, term1)
        Arblib.add!(S1, S1, abs_term)
        Arblib.abs!(abs_term, term2)
        Arblib.add!(S2, S2, abs_term)
    end

    C_R_U_1, C_R_U_2 = C_R_U_12(n, a, b, z₁)

    # Note that Γ'(a) / Γ(a) is exactly the digamma function
    R = (1 + abs(digamma(a) / log(z₁))) * C_R_U(n, a, b, z₁) + C_R_U_1 + C_R_U_2

    return S1 + S2 / abs(log(z₁)) + R * abs(z₁)^-n
end
