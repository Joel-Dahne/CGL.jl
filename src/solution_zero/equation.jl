export cgl_equation_real, cgl_equation_real_alt, cgl_equation_real_taylor

"""
    ivp_zero_complex(μ, κ, λ)

The initial value problem given by
```
0 = (1 - im * ϵ) * (d2(Q) + (d - 1) / ξ * d(Q)) +
    im * κ * ξ * d(Q) + im * κ / σ * Q - ω * Q + (1 + im * δ) * abs(Q)^2σ * Q
```
with
```
Q(0) = μ
d(Q)(0) = 0
```
We here use `d(Q)` and `d2(Q)` to denote the first and second
derivative of `Q` respectively.
"""
function ivp_zero_complex(μ, κ, λ) end

"""
    ivp_zero_real(μ, κ, λ)

The initial value problem given by
```
TODO
```
with
```
a(0) = μ
d(a)(0) = 0
b(0) = 0
d(b)(0) = 0
```
We here use `d(a)` and `d2(a)` to denote the first and second
derivative of `a` respectively.
"""
function ivp_zero_real(μ, κ, λ) end

"""
    ivp_zero_real_system(μ, κ, λ)

The initial value problem given by
```
d(a) = α
d(b) = β
d(α) = F1(a, b, α, β, ξ) - ϵ * F2(a, b, α, β, ξ)
d(β) = ϵ * F1(a, b, α, β, ξ) + F2(a, b, α, β, ξ)
```
with
```
a(0) = μ
b(0) = 0
α(0) = 0
β(0) = 0
```
We here use `d(a)` to denote the derivative of `a`. Here `F1` and `F2`
are given by
```
F1(a, b, α, β, ξ) = -(d - 1) / ξ * (α + ϵ * β) +
    κ * ξ * β +
    κ / σ * b +
    ω * a -
    (a^2 + b^2)^σ * a +
    δ * (a^2 + b^2)^σ * b
```
and
```
F2(a, b, α, β, ξ) = -(d - 1) / ξ * (α - ϵ * β) -
    κ * ξ * α -
    κ / σ * a +
    ω * b -
    (a^2 + b^2)^σ * b -
    δ * (a^2 + b^2)^σ * a
```
"""
function ivp_zero_real_system(μ, κ, λ) end

"""
    ivp_zero_real_system_autonomus(μ, κ, λ)

The initial value problem given by
```
d(a) = α
d(b) = β
d(α) = F1(a, b, α, β, ξ) - ϵ * F2(a, b, α, β, ξ)
d(β) = ϵ * F1(a, b, α, β, ξ) + F2(a, b, α, β, ξ)
d(ξ) = 1
```
with
```
a(0) = μ
b(0) = 0
α(0) = 0
β(0) = 0
ξ(0) = 0
```
We here use `d(a)` to denote the derivative of `a`. `F1` and `F2` are
as in [`ivp_zero_real_system`](@ref).
"""
function ivp_zero_real_system_autonomus(μ, κ, λ) end

"""
    cgl_equation_real(u, κ, ξ, λ)

Evaluate the right hand side of [`ivp_zero_real_system`](@ref) at the
point `u = [a, b, α, β]` and time `ξ`.

For `ξ = 0` the first term for both `F1` and `F2` have a division by
zero. For this term to be finite we in this case need `α = β = 0`, if
that is the case we set the term to zero.
- **TODO:** Is fixing the term to be zero the correct thing to do?
"""
function cgl_equation_real(u, κ, ξ, λ::CGLParams)
    a, b, α, β = u
    d, ω, σ, ϵ, δ = λ.d, λ.ω, λ.σ, λ.ϵ, λ.δ

    a2b2σ = (a^2 + b^2)^σ

    # TODO: Should we keep fastmath? It makes it faster and probably
    # improves accuracy due to fma.
    @fastmath F1 = κ * ξ * β + κ / σ * b + ω * a - a2b2σ * a + δ * a2b2σ * b
    @fastmath F2 = -κ * ξ * α - κ / σ * a + ω * b - a2b2σ * b - δ * a2b2σ * a

    if !isone(d) && !(iszero(ξ) && iszero(α) && iszero(β))
        F1 -= (d - 1) / ξ * (α + ϵ * β)
        F2 -= (d - 1) / ξ * (β - ϵ * α)
    end

    return SVector(α, β, F1 - ϵ * F2, ϵ * F1 + F2)
end

"""
    cgl_equation_real_alt(u, ξ, (κ, λ))

Like [`cgl_equation_real`](@ref) but with an interface that works for
[`ODEProblem`](@ref).
"""
cgl_equation_real_alt(u, (κ, λ), ξ) = cgl_equation_real(u, κ, ξ, λ)

"""
    cgl_equation_real_taylor(((a0, a1), (b0, b1)), κ, ξ₀, λ; degree = 5)

Compute the expansion of `a` and `b` in [`cgl_equation_real`](@ref)
centered at the point `ξ = ξ₀` with the first two coefficients in the
expansions for `a` and `b` given by `a0, a1` and `b0, b1`
respectively.
"""
function cgl_equation_real_taylor(
    u0::AbstractVector{NTuple{2,Arb}},
    κ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree::Integer = 5,
)
    d, ω, σ, ϵ, δ = λ.d, λ.ω, λ.σ, λ.ϵ, λ.δ

    a = ArbSeries(u0[1]; degree)
    b = ArbSeries(u0[2]; degree)

    if iszero(ξ₀) && !isone(d) && !iszero(a[1]) && !iszero(b[1])
        return SVector(indeterminate(a), indeterminate(b))
    end

    for n = 0:degree-2
        a2b2σ = (ArbSeries(a, degree = n + 1)^2 + ArbSeries(b, degree = n + 1)^2)^σ
        u1 = a2b2σ * a
        u2 = a2b2σ * b

        if iszero(ξ₀)
            F1 = κ * n * b[n] + κ / σ * b[n] + ω * a[n] - u1[n] + δ * u2[n]
            F2 = -κ * n * a[n] - κ / σ * a[n] + ω * b[n] - u2[n] - δ * u1[n]

            a[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + d))
            b[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + d))
        else

            F1 =
                κ * (ξ₀ * (n + 1) * b[n+1] + n * b[n]) + κ / σ * b[n] + ω * a[n] - u1[n] +
                δ * u2[n]
            F2 =
                -κ * (ξ₀ * (n + 1) * a[n+1] + n * a[n]) - κ / σ * a[n] + ω * b[n] - u2[n] -
                δ * u1[n]

            if !isone(d)
                v1 = Arblib.derivative(a) / ArbSeries((ξ₀, 1), degree = n) # a' / ξ
                v2 = Arblib.derivative(b) / ArbSeries((ξ₀, 1), degree = n) # b' / ξ
                F1 += (1 - d) * (v1[n] + ϵ * v2[n])
                F2 += (1 - d) * (v2[n] - ϵ * v1[n])
            end

            a[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + 1))
            b[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + 1))
        end
    end

    return SVector(a, b)
end

"""
    cgl_equation_real_dμ_taylor(((a_dμ0, a_dμ1), (b_dμ0, b_dμ1)), κ, ξ₀, λ; degree = 5)

Compute the expansion of `a_dμ` and `b_dμ` in
[`cgl_equation_real`](@ref) centered at the point `ξ = ξ₀` with the
first two coefficients in the expansions for `a_dμ` and `b_dμ` given
by `a0_dμ, a1_dμ` and `b0_dμ, b1_dμ` respectively.
"""
function cgl_equation_real_dμ_taylor(
    u0_dμ::AbstractVector{NTuple{2,Arb}},
    a::ArbSeries,
    b::ArbSeries,
    κ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree::Integer = 5,
)
    @assert Arblib.degree(a) == Arblib.degree(b) == degree

    d, ω, σ, ϵ, δ = λ.d, λ.ω, λ.σ, λ.ϵ, λ.δ

    a_dμ = ArbSeries(u0_dμ[1]; degree)
    b_dμ = ArbSeries(u0_dμ[2]; degree)

    if iszero(ξ₀) && !isone(d) && !iszero(a_dμ[1]) && !iszero(b_dμ[1])
        return SVector(indeterminate(a_dμ), indeterminate(b_dμ))
    end

    a2b2σ = (a^2 + b^2)^σ
    d_a2b2σ_aa = 2σ * (a^2 + b^2)^(σ - 1) * a^2
    d_a2b2σ_ab = 2σ * (a^2 + b^2)^(σ - 1) * a * b
    d_a2b2σ_bb = 2σ * (a^2 + b^2)^(σ - 1) * b^2

    for n = 0:degree-2
        u1 =
            a2b2σ * ArbSeries(a_dμ, degree = n + 1) +
            d_a2b2σ_aa * ArbSeries(a_dμ, degree = n + 1) +
            d_a2b2σ_ab * ArbSeries(b_dμ, degree = n + 1)
        u2 =
            a2b2σ * ArbSeries(b_dμ, degree = n + 1) +
            d_a2b2σ_ab * ArbSeries(a_dμ, degree = n + 1) +
            d_a2b2σ_bb * ArbSeries(b_dμ, degree = n + 1)

        if iszero(ξ₀)
            F1 = κ * n * b_dμ[n] + κ / σ * b_dμ[n] + ω * a_dμ[n] - u1[n] + δ * u2[n]
            F2 = -κ * n * a_dμ[n] - κ / σ * a_dμ[n] + ω * b_dμ[n] - u2[n] - δ * u1[n]

            a_dμ[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + d))
            b_dμ[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + d))
        else
            error("not implemented")
        end
    end

    return SVector(a_dμ, b_dμ)
end

"""
    cgl_equation_real_dκ_taylor(((a_dμ0, a_dμ1), (b_dμ0, b_dμ1)), κ, ξ₀, λ; degree = 5)

Compute the expansion of `a_dκ` and `b_dκ` in [`equation_real`](@ref)
centered at the point `ξ = ξ₀` with the first two coefficients in the
expansions for `a_dκ` and `b_dκ` given by `a0_dκ, a1_dκ` and `b0_dκ,
b1_dκ` respectively.
"""
function cgl_equation_real_dκ_taylor(
    u0_dκ::AbstractVector{NTuple{2,Arb}},
    a::ArbSeries,
    b::ArbSeries,
    κ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree::Integer = 5,
)
    @assert Arblib.degree(a) == Arblib.degree(b) == degree

    d, ω, σ, ϵ, δ = λ.d, λ.ω, λ.σ, λ.ϵ, λ.δ

    a_dκ = ArbSeries(u0_dκ[1]; degree)
    b_dκ = ArbSeries(u0_dκ[2]; degree)

    if iszero(ξ₀) && !isone(d) && !iszero(a_dκ[1]) && !iszero(b_dκ[1])
        return SVector(indeterminate(a_dκ), indeterminate(b_dκ))
    end

    a2b2σ = (a^2 + b^2)^σ
    d_a2b2σ_aa = 2σ * (a^2 + b^2)^(σ - 1) * a^2
    d_a2b2σ_ab = 2σ * (a^2 + b^2)^(σ - 1) * a * b
    d_a2b2σ_bb = 2σ * (a^2 + b^2)^(σ - 1) * b^2

    for n = 0:degree-2
        u1 =
            a2b2σ * ArbSeries(a_dκ, degree = n + 1) +
            d_a2b2σ_aa * ArbSeries(a_dκ, degree = n + 1) +
            d_a2b2σ_ab * ArbSeries(b_dκ, degree = n + 1)
        u2 =
            a2b2σ * ArbSeries(b_dκ, degree = n + 1) +
            d_a2b2σ_ab * ArbSeries(a_dκ, degree = n + 1) +
            d_a2b2σ_bb * ArbSeries(b_dκ, degree = n + 1)

        if iszero(ξ₀)
            F1 =
                κ * n * b_dκ[n] + n * b[n] + κ / σ * b_dκ[n] + 1 / σ * b[n] + ω * a_dκ[n] -
                u1[n] + δ * u2[n]
            F2 =
                -κ * n * a_dκ[n] - n * a[n] - κ / σ * a_dκ[n] - 1 / σ * a[n] + ω * b_dκ[n] -
                u2[n] - δ * u1[n]

            a_dκ[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + d))
            b_dκ[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + d))
        else
            error("not implemented")
        end
    end

    return SVector(a_dκ, b_dκ)
end
