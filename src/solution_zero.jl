export gl_equation_real

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
function gl_equation_real(u, (p, κ)::Tuple{GLParams{T},T}, ξ) where {T}
    a, b, α, β = u

    d, ω, σ, ϵ, δ = p.d, p.ω, p.σ, p.ϵ, p.δ

    a2b2σ = (a^2 + b^2)^σ

    F1 = κ * ξ * β + κ / σ * b + ω * a - a2b2σ * a + δ * a2b2σ * b
    F2 = -κ * ξ * α - κ / σ * a + ω * b - a2b2σ * b - δ * a2b2σ * a

    if !(iszero(ξ) && iszero(α) && iszero(β))
        F1 -= (d - 1) / ξ * (α + ϵ * β)
        F2 -= (d - 1) / ξ * (β - ϵ * α)
    end

    return SVector(α, β, F1 - ϵ * F2, ϵ * F1 + F2)
end
