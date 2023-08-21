"""
    fpp_infinity_complex(γ, κ, λ)

The fixed point problem given by
```
Q(ξ) = T(Q)(ξ)
```
where `T(Q)` is the operator
```
T(Q) = γ * P(ξ) - ∫_ξ₁^∞ (1 + im * δ) * K(ξ, η) * abs(Q(η))^2σ * Q(η) dη
```
"""
function fpp_infinity_complex(γ, κ, λ) end
