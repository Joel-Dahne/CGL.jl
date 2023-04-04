export gl_equation_real, gl_taylor_expansion_real

"""
    gl_equation_real(u, p, ξ)

In real variables the equation can be written as a system of first
order ODE:s in the four variables `a, b, α, β`. The equation is then
given by
```
deriv(a) = α
deriv(b) = β
deriv(α) = F1(a, b, α, β, ξ) - ϵ * F2(a, b, α, β, ξ)
deriv(β) = ϵ * F1(a, b, α, β, ξ) + F2(a, b, α, β, ξ)
```
where we use `deriv(a)` to denote the derivative of `a`. Here
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
For `ξ = 0` the first term for both `F1` and `F2` have a division by
zero. For this term to be finite we in this case need `α = β = 0`, if
that is the case we set the term to zero.
- **TODO:** Is fixing the term to be zero the correct thing to do?
"""
function gl_equation_real(u, (p, κ)::Tuple{AbstractGLParams{T},T}, ξ) where {T}
    a, b, α, β = u

    d, ω, σ, ϵ, δ = p.d, p.ω, p.σ, p.ϵ, p.δ

    a2b2σ = (a^2 + b^2)^σ

    F1 = κ * ξ * β + κ / σ * b + ω * a - a2b2σ * a + δ * a2b2σ * b
    F2 = -κ * ξ * α - κ / σ * a + ω * b - a2b2σ * b - δ * a2b2σ * a

    if !isone(d) && !(iszero(ξ) && iszero(α) && iszero(β))
        F1 -= (d - 1) / ξ * (α + ϵ * β)
        F2 -= (d - 1) / ξ * (β - ϵ * α)
    end

    return SVector(α, β, F1 - ϵ * F2, ϵ * F1 + F2)
end

"""
    gl_taylor_expansion_real((a0, b0, a1, b1), ξ0, degree , (p, κ))

Compute the expansion of `a` and `b` centered at the point `ξ = ξ0`
with the first two coefficients in the expansions for `a` and `b`
given by `a0, a1` and `b0, b1` respectively.
"""
function gl_taylor_expansion_real(
    u0::AbstractVector{Arb},
    ξ0::Arb,
    (p, κ)::Tuple{AbstractGLParams{Arb},Arb};
    degree::Integer = 5,
)
    a0, b0, a1, b1 = u0

    d, ω, σ, ϵ, δ = p.d, p.ω, p.σ, p.ϵ, p.δ

    a = ArbSeries((a0, a1); degree)
    b = ArbSeries((b0, b1); degree)

    if iszero(ξ0) && !isone(d) && !iszero(a1) && !iszero(b1)
        return indeterminate(a), indeterminate(b)
    end

    for n = 0:degree-2
        a2b2σ = (a^2 + b^2)^σ
        u1 = a2b2σ * a
        u2 = a2b2σ * b

        if iszero(ξ0)
            F1 = κ * n * b[n] + κ / σ * b[n] + ω * a[n] - u1[n] + δ * u2[n]
            F2 = -κ * n * a[n] - κ / σ * a[n] + ω * b[n] - u2[n] - δ * u1[n]

            a[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + d))
            b[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + d))
        else
            F1 =
                κ * (ξ0 * (n + 1) * b[n+1] + n * b[n]) + κ / σ * b[n] + ω * a[n] - u1[n] +
                δ * u2[n]
            F2 =
                -κ * (ξ0 * (n + 1) * a[n+1] + n * a[n]) - κ / σ * a[n] + ω * b[n] - u2[n] -
                δ * u1[n]

            if !isone(d)
                v1 = Arblib.derivative(a) / ArbSeries((ξ0, 1), degree = n) # a' / ξ
                v2 = Arblib.derivative(b) / ArbSeries((ξ0, 1), degree = n) # b' / ξ
                F1 += (1 - d) * (v1[n] + ϵ * v2[n])
                F2 += (1 - d) * (v2[n] - ϵ * v1[n])
            end

            a[n+2] = (F1 - ϵ * F2) / ((n + 2) * (n + 1))
            b[n+2] = (ϵ * F1 + F2) / ((n + 2) * (n + 1))
        end
    end

    return a, b
end
