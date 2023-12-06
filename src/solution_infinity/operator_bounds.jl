"""
    C_T1(v::Arb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)

Return `C₁` and `C₂` such that
```
norm(T(u, (γ, κ, p)), v) <= C₁ * abs(γ) * ξ₁^(-v) + C₂ * ξ₁^(-2 + 2σ * v) * norm(u, v)^(2σ + 1)
```
Here `norm(u, v)` is given by ``\\sup_{ξ ≥ ξ₁} |ξ|^{1 / σ - v}|u(ξ)|``
and `T` is the operator
```
T(u, (γ, κ, p))(ξ) = γ * P(ξ, (p, κ)) - ∫ (1 + im * δ) * K(ξ, η, (p, κ)) * abs(u(η))^2σ * u(η) dη
```
with the integration taken from `ξ₁` to `∞`. We require that
```
(2σ + 1)v < 2 + 2 / σ - d
```

We have
```
abs(ξ)^(1 / σ - v) * abs(T(u)) <= abs(γ) * abs(ξ)^(1 / σ - v) * abs(P(ξ)) +
    abs(1 + im * δ) * abs(ξ)^(1 / σ - v) * ∫ abs(K(ξ, η, (p, κ))) * abs(u(η))^(2σ + 1) dη
```

For the first term we have
```
abs(γ) * abs(ξ)^(1 / σ - v) * abs(P(ξ)) <= CP * abs(γ) * ξ^(-v) <= CP * abs(γ) * ξ₁^(-v)
```
with `CP` as given by [`C_P`](@ref). So we can take `C₁ = CP`.

For the second term we use that
```
abs(u(η))^(2σ + 1) <= norm(abs(u(η))^(2σ + 1), v)^(2σ + 1) * abs(η)^((2σ + 1) * (v - 1 / σ))
```
to get the bound
```
abs(1 + im * δ) * abs(ξ)^(1 / σ - v) * ∫ abs(K(ξ, η, (p, κ))) * abs(u(η))^(2σ + 1) dη <=
    abs(1 + im * δ) * norm(abs(u), v)^(2σ + 1) * abs(ξ)^(1 / σ - v) *
        ∫ abs(K(ξ, η, (p, κ))) * abs(η)^((2σ + 1) * (v - 1 / σ)) dη
```
Splitting the integral at `ξ` and letting `CK` be the bound from
[`C_K`](@ref) we get
```
abs(ξ)^(1 / σ - v) * ∫ abs(K(ξ, η, (p, κ))) * abs(η)^((2σ + 1) * (v - 1 / σ)) dη <=
    CK * (
        ξ^(-v) * ∫_ξ₁^ξ η^(-3 + (2σ + 1) * v) dη +
        ξ^(2 / σ - d - v) * ∫_ξ^∞ η^(-3 + (2σ + 1) * v + d - 2 / σ) dη
)
```
Since `(2σ + 1)v < 2 + 2 / σ - d` we have that the integrals converge
and are bounded by
```
ξ^(-v) * ∫_ξ₁^ξ η^(-3 + (2σ + 1) * v) dη <= 2ξ₁^(-2 + 2σ * v) / abs(-2 + (2σ + 1) * v)

ξ^(2 / σ - d - v) * ∫_ξ^∞ η^(-3 + (2σ + 1) * v + d - 2 / σ) dη <=
    ξ₁^(-2 + 2σ * v) / abs(-2 + (2σ + 1) * v + d - 2 / σ)
```
Hence the second term is bounded by
```
abs(1 + im * δ) * CK * norm(abs(u), v)^(2σ + 1) * ξ^(-2 + 2σ * v) * (
    2 / abs(-2 + (2σ + 1) * v) + 1 / abs(-2 + (2σ + 1) * v + d - 2 / σ)
)
```
and we can take
```
C₂ = abs(1 + im * δ) * CK * (
    2 / abs(-2 + (2σ + 1) * v) + 1 / abs(-2 + (2σ + 1) * v + d - 2 / σ)
)
```

"""
function C_T1(v::Arb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)
    d, ω, σ, ϵ, δ = p.d, p.ω, p.σ, p.ϵ, p.δ

    (2σ + 1)v < 2 + 2 / σ - d || throw(ArgumentError("assuming (2σ + 1)v < 2 + 2 / σ - d"))


    # TODO: This might follow from the other requirements?
    @assert -2 + (2σ + 1) * v < 0
    @assert -2 + (2σ + 1) * v + d - 2 / σ < 0

    C₁ = C_P(κ, p, ξ₁)

    C₂ =
        abs(1 + im * δ) *
        C_K(κ, p, ξ₁) *
        (1 / abs(-2 + (2σ + 1) * v) + 1 / abs(-2 + (2σ + 1) * v + d - 2 / σ))

    return C₁, C₂
end

"""
    C_T2(v::Arb, γ::Acb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)

Return `C` such that
```
norm(T(u₁, (γ, κ, p)) - T(u₂, (γ, κ, p)), v) <=
    C * ξ₁^(-2 + 2σ * v) * norm(u₁ - u₂, v) * (norm(u₁, v)^2σ + norm(u₂, v)^2σ)
```
Here `norm(u, v)` and T(u, (γ, κ, p)) are defined as in
[`C_T1`](@ref).

# Helper inequality
We make use of the inequality
```
abs(abs(z₁)^2σ * z₁ - abs(z₂)^2σ * z₂) <= M * abs(z₁ - z₂) * (abs(z₁)^2σ + abs(z₂)^2σ)
```
that holds for all complex `z₁` and `z₂` with properly chosen `M`.

# Finding `C`
We have
```
ξ^(1 / σ - v) * abs(T(u₁) - T(u₂)) <= abs(1 + im * δ) * ξ^(1 / σ - v) *
    ∫abs(K(ξ, η)) * abs(abs(u₁(η))^2σ * u₁(η) - abs(u₂(η))^2σ * u₂(η)) dη
```
From the helper inequality we get
```
∫abs(K(ξ, η)) * abs(abs(u₁(η))^2σ * u₁(η) - abs(u₂(η))^2σ * u₂(η)) dη <=
    M * ∫abs(K(ξ, η)) * abs(u₁(η) - u₂(η)) * abs(abs(u₁(η))^2σ - abs(u₂(η))^2σ) dη <=
    M * norm(u₁ - u₂, v) * (norm(u₁, v)^2σ + norm(u₂, v)^2σ) * ∫abs(K(ξ, η)) * η^((2σ + 1) * (v - 1 / σ)) dη
```
Hence
```
ξ^(1 / σ - v) * abs(T(u₁) - T(u₂)) <= M * abs(1 + im * δ) * ξ^(1 / σ - v) * norm(u₁ - u₂, v) *
    (norm(u₁, v)^2σ + norm(u₂, v)^2σ) * ∫abs(K(ξ, η)) * η^((2σ + 1) * (v - 1 / σ)) dη
```
With `C₂` as in [`C_T1`](@ref) we get that this is bounded by
```
M * C₂ * ξ₁^(-2 + 2σ * v) * norm(u₁ - u₂, v) * (norm(u₁, v)^2σ + norm(u₂, v)^2σ)
```
and we see that we can take `C = M * C₂`.
"""
function C_T2(v::Arb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)
    M = if isone(p.σ)
        sqrt(Arb(2)) / (4 - 2sqrt(Arb(2))) + 1
    else
        (; σ) = p

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

        max(M1, M2, M3, M4)
    end

    return M * C_T1(v, κ, p, ξ₁)[2]
end
