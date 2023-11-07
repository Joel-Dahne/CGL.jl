export gl_equation_real_system,
    gl_equation_real_system_ode,
    gl_equation_real_taylor_expansion,
    gl_equation_real_system_autonomus_taylor_expansion,
    gl_equation_real_system_autonomus_taylor_expansion_simple

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
    gl_equation_real_system(u, κ, ξ, λ)

Evaluate the right hand side of [`ivp_zero_real_system`](@ref) at the
point `u = [a, b, α, β]` and time `ξ`.

For `ξ = 0` the first term for both `F1` and `F2` have a division by
zero. For this term to be finite we in this case need `α = β = 0`, if
that is the case we set the term to zero.
- **TODO:** Is fixing the term to be zero the correct thing to do?
"""
function gl_equation_real_system(u, κ, ξ, λ::AbstractGLParams)
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
    gl_equation_real_system_ode(u, ξ, (κ, λ))

Like [`gl_equation_real_system`](@ref) but with an interface that
works for [`ODEProblem`](@ref).
"""
gl_equation_real_system_ode(u, (κ, λ), ξ) = gl_equation_real_system(u, κ, ξ, λ)

"""
    gl_equation_real_taylor_expansion(((a0, a1), (b0, b1)), κ, ξ₀, λ; degree = 5)

Compute the expansion of `a` and `b` in [`equation_real`](@ref)
centered at the point `ξ = ξ₀` with the first two coefficients in the
expansions for `a` and `b` given by `a0, a1` and `b0, b1`
respectively.
"""
function gl_equation_real_taylor_expansion(
    u0::AbstractVector{NTuple{2,Arb}},
    κ::Arb,
    ξ₀::Arb,
    λ::AbstractGLParams{Arb};
    degree::Integer = 5,
)
    d, ω, σ, ϵ, δ = λ.d, λ.ω, λ.σ, λ.ϵ, λ.δ

    a = ArbSeries(u0[1]; degree)
    b = ArbSeries(u0[2]; degree)

    if iszero(ξ₀) && !isone(d) && !iszero(a[1]) && !iszero(a[1])
        return [indeterminate(a), indeterminate(b)]
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
    gl_equation_real_system_autonomus_taylor_expansion((a0, b0, α0, β0, ξ0), κ, λ; degree)

Compute the expansion of `a, b, α, β, ξ` in
[`equation_real_system_autonomus`](@ref), centered at zero.
- **TODO:** Where exactly is the expansion computed?
"""
function gl_equation_real_system_autonomus_taylor_expansion(
    u0::AbstractVector{Arb},
    κ::Arb,
    λ::AbstractGLParams{Arb};
    degree::Integer = 5,
)
    a0, b0, α0, β0, ξ0 = u0

    d, ω, σ, ϵ, δ = λ.d, λ.ω, λ.σ, λ.ϵ, λ.δ

    σ_inv = inv(σ)

    a = ArbSeries(a0; degree)
    b = ArbSeries(b0; degree)
    α = ArbSeries(α0; degree)
    β = ArbSeries(β0; degree)
    ξ = ArbSeries((ξ0, 1); degree)

    if iszero(ξ0) && !isone(d) && !iszero(α0) && !iszero(β0)
        return [
            indeterminate(ξ),
            indeterminate(a),
            indeterminate(b),
            indeterminate(α),
            indeterminate(β),
        ]
    end

    # Holds nth coefficient of u1 = (a^2 + b^2)^σ * a and u2 = (a^2 +
    # b^2)^σ * b
    u1n, u2n = zero(Arb), zero(Arb)

    # Used internally by _a2b2σ!
    a2b2 = ArbRefVector(degree)
    a2b2_diff = ArbRefVector(degree - 1)
    a2b2_inv = ArbRefVector(degree - 1)
    a2b2_div = ArbRefVector(degree - 1)
    σ_log_a2b2 = ArbRefVector(degree)
    a2b2σ = ArbRefVector(degree)

    F1, F2 = zero(Arb), zero(Arb)

    if !isone(d)
        v1 = zero(α)
        v2 = zero(β)
        tmp = zero(Arb)
    end

    for n = 0:degree-1
        _u1nu2n!(
            u1n,
            u2n,
            a2b2,
            a2b2_diff,
            a2b2_inv,
            a2b2_div,
            σ_log_a2b2,
            a2b2σ,
            a,
            b,
            σ,
            n,
        )

        # Compute
        # F1 = κ * ((ξ * β)[n] + b[n] / σ) + ω * a[n] - u1n + δ * u2n
        # F2 = -κ * ((ξ * α)[n] + a[n] / σ) + ω * b[n] - u2n - δ * u1n
        # Compute (ξ*β)[n] and (ξ*α)[n]
        if iszero(n)
            Arblib.mul!(F1, ξ0, Arblib.ref(β, n))
            Arblib.mul!(F2, ξ0, Arblib.ref(α, n))
        else
            Arblib.fma!(F1, ξ0, Arblib.ref(β, n), Arblib.ref(β, n - 1))
            Arblib.fma!(F2, ξ0, Arblib.ref(α, n), Arblib.ref(α, n - 1))
        end

        # IMPROVE: Consider doing a division instead of
        # multiplication, it might give slightly better enclosure but
        # is slower
        # IMPROVE: It seems like it might be slightly better to
        # multiply the terms by κ independently. It gives marginally
        # better enclosures if κ is exact at least.
        Arblib.addmul!(F1, σ_inv, Arblib.ref(b, n))
        Arblib.addmul!(F2, σ_inv, Arblib.ref(a, n))

        Arblib.mul!(F1, F1, κ)
        Arblib.mul!(F2, F2, κ)
        Arblib.neg!(F2, F2)

        Arblib.addmul!(F1, ω, Arblib.ref(a, n))
        Arblib.addmul!(F2, ω, Arblib.ref(b, n))

        Arblib.sub!(F1, F1, u1n)
        Arblib.sub!(F2, F2, u2n)

        Arblib.addmul!(F1, δ, u2n)
        Arblib.submul!(F2, δ, u1n)

        if iszero(ξ0)
            # Setting it like this ensures that the degree gets
            # updated correctly
            α[n+1] = F1
            β[n+1] = F2

            Arblib.submul!(Arblib.ref(α, n + 1), ϵ, F2)
            Arblib.addmul!(Arblib.ref(β, n + 1), ϵ, F1)

            Arblib.div!(Arblib.ref(α, n + 1), Arblib.ref(α, n + 1), n + d)
            Arblib.div!(Arblib.ref(β, n + 1), Arblib.ref(β, n + 1), n + d)
        else
            if !isone(d)
                # IMPROVE: This currently computes the full series every time
                Arblib.div_series!(v1, α, ξ, n + 1)
                Arblib.div_series!(v2, β, ξ, n + 1)

                # F1 += (1 - d) * (v1[n] + ϵ * v2[n])
                Arblib.fma!(tmp, ϵ, Arblib.ref(v2, n), Arblib.ref(v1, n))
                Arblib.addmul!(F1, tmp, 1 - d)

                # F2 += (1 - d) * (v2[n] - ϵ * v1[n])
                Arblib.set!(tmp, Arblib.ref(v2, n))
                Arblib.submul!(tmp, ϵ, Arblib.ref(v1, n))
                Arblib.addmul!(F2, tmp, 1 - d)
            end

            # Setting it like this ensures that the degree gets
            # updated correctly
            α[n+1] = F1
            β[n+1] = F2

            Arblib.submul!(Arblib.ref(α, n + 1), ϵ, F2)
            Arblib.addmul!(Arblib.ref(β, n + 1), ϵ, F1)

            Arblib.div!(Arblib.ref(α, n + 1), Arblib.ref(α, n + 1), n + 1)
            Arblib.div!(Arblib.ref(β, n + 1), Arblib.ref(β, n + 1), n + 1)
        end

        a[n+1] = Arblib.ref(α, n)
        b[n+1] = Arblib.ref(β, n)

        Arblib.div!(Arblib.ref(a, n + 1), Arblib.ref(a, n + 1), n + 1)
        Arblib.div!(Arblib.ref(b, n + 1), Arblib.ref(b, n + 1), n + 1)
    end

    return SVector(a, b, α, β, ξ)
end


"""
    gl_equation_real_system_autonomus_taylor_expansion_simple((a0, b0, α0, β0, ξ0), κ, λ; degree)

Same as [`gl_equation_real_system_autonomus_taylor_expansion`](@ref)
but less optimized, mainly used for testing.
"""
function gl_equation_real_system_autonomus_taylor_expansion_simple(
    u0::AbstractVector{Arb},
    κ::Arb,
    λ::AbstractGLParams{Arb};
    degree::Integer = 5,
)
    a0, b0, α0, β0, ξ0 = u0

    d, ω, σ, ϵ, δ = λ.d, λ.ω, λ.σ, λ.ϵ, λ.δ

    a = ArbSeries(a0, ; degree)
    b = ArbSeries(b0; degree)
    α = ArbSeries(α0; degree)
    β = ArbSeries(β0; degree)
    ξ = ArbSeries((ξ0, 1); degree)

    if iszero(ξ0) && !isone(d) && !iszero(α0) && !iszero(β0)
        return [
            indeterminate(ξ),
            indeterminate(a),
            indeterminate(b),
            indeterminate(α),
            indeterminate(β),
        ]
    end

    for n = 0:degree-1
        a2b2σ = (ArbSeries(a, degree = n)^2 + ArbSeries(b, degree = n)^2)^σ

        u1 = a2b2σ * a
        u2 = a2b2σ * b

        if iszero(ξ0)
            F1 = κ * (ξ*β)[n] + κ / σ * b[n] + ω * a[n] - u1[n] + δ * u2[n]
            F2 = -κ * (ξ*α)[n] - κ / σ * a[n] + ω * b[n] - u2[n] - δ * u1[n]

            a[n+1] = α[n] / (n + 1)
            b[n+1] = β[n] / (n + 1)
            α[n+1] = (F1 - ϵ * F2) / (n + d)
            β[n+1] = (ϵ * F1 + F2) / (n + d)
        else
            # Compute (ξ*β)[n] and (ξ*α)[n]
            if iszero(n)
                ξβn = ξ[0] * β[n]
                ξαn = ξ[0] * α[n]
            else
                ξβn = ξ[0] * β[n] + ξ[1] * β[n-1]
                ξαn = ξ[0] * α[n] + ξ[1] * α[n-1]
            end

            F1 = κ * ξβn + κ / σ * b[n] + ω * a[n] - u1[n] + δ * u2[n]
            F2 = -κ * ξαn - κ / σ * a[n] + ω * b[n] - u2[n] - δ * u1[n]

            if !isone(d)
                v1 = α / ξ
                v2 = β / ξ
                F1 += (1 - d) * (v1[n] + ϵ * v2[n])
                F2 += (1 - d) * (v2[n] - ϵ * v1[n])
            end

            a[n+1] = α[n] / (n + 1)
            b[n+1] = β[n] / (n + 1)
            α[n+1] = (F1 - ϵ * F2) / (n + 1)
            β[n+1] = (ϵ * F1 + F2) / (n + 1)
        end
    end

    return SVector(a, b, α, β, ξ)
end
