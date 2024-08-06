"""
    fpp_infinity_complex(γ, κ, λ)

The fixed point problem given by
```
Q(ξ) = T(Q)(ξ)
```
where `T(u)` is the operator
```
T(u) = γ * P(ξ) +
    P(ξ) * ∫_ξ₁^ξ  J_E(ξ) * abs(u(η))^2σ * u(η) dη +
    E(ξ) * ∫_ξ^∞  J_P(ξ) * abs(u(η))^2σ * u(η) dη
```
"""
function fpp_infinity_complex(γ, κ, λ) end
