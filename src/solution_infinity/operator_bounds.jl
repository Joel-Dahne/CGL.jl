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
        (2 / abs(-2 + (2σ + 1) * v) + 1 / abs(-2 + (2σ + 1) * v + d - 2 / σ))

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

**FIXME:** Currently we use a hard coded value for `M`.

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
    # FIXME
    M = if isone(p.σ)
        (4 - sqrt(Arb(2))) / (4 - 2sqrt(Arb(2))) - 1
    elseif p.σ ≈ 2.3
        Arb(3.383) # FIXME
    else
        error("no approximate M for given σ")
    end

    return M * C_T1(v, κ, p, ξ₁)[2]
end

"""
    C_fix_point(r₁::Arb, v::Arb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)

Return `C₁` and `C₂` such that the fixed point `u(ξ)` of the operator
`T` satisfies
```
u(ξ) = γ * P(ξ) + P(ξ) * f(ξ) + E(ξ) * g(ξ)
```
where
```
abs(f(ξ)) <= C₁ * (ξ₁^(2σ * v + v - 2) - ξ^(2σ * v + v - 2))
abs(g(ξ)) <= C₂ * exp(-real(c) * ξ^2) * ξ^(-2 / σ + 2σ * v + v + d - 2)
```
for `ξ >= ξ₁` and `abs(γ) <= r₁`. In particular for `ξ = ξ₁` we have
```
u(ξ) = γ * P(ξ) + E(ξ) * g(ξ)
```
Here
```
real(c) = κ * ϵ / (1 + ϵ^2)
```

# Splitting `T`
Since `u` is a fix point of `T` we have
```
u(ξ) = γ * P(ξ) - (1 + im * δ) * ∫_ξ₁^∞ K(ξ, η) * abs(u(η))^2σ * u(η) dη
```
If we let
```
f(ξ) = (1 + im * δ) / (1 - im * ϵ) * ∫_ξ₁^ξ E(η) / W(η) * abs(u(η))^2σ * u(η) dη

g(ξ) = (1 + im * δ) / (1 - im * ϵ) * ∫_ξ^∞ P(η) / W(η) * abs(u(η))^2σ * u(η) dη
```
we can write the above as
```
u(ξ) = γ * P(ξ) + P(ξ) * f(ξ) + E(ξ) * g(ξ)
```

# Bounding `abs(f(ξ))`
Focusing on the integral we have
```
abs(∫_ξ₁^ξ E(η) / W(η) * abs(u(η))^2σ * u(η) dη) <=
    ∫_ξ₁^ξ abs(E(η)) / abs(W(η)) * abs(u(η))^(2σ + 1) dη
```
Using that
```
abs(E(η)) <= CE * exp(real(c) * η^2) * η^(-d + 1 / σ)
abs(W(η)) =  abs(2c) * exp(-imag(b - a) * π) * abs(c^-b) * η^(1 - d) * exp(real(c) * η^2)
abs(u(η)) <= norm(u, v) * η^(v - 1 / σ)
```
we get
```
abs(∫_ξ₁^ξ E(η) / W(η) * abs(u(η))^2σ * u(η) dη) <=
    CE / abs(2c) * exp(imag(b - a) * π) * abs(c^b) * norm(u, v)^(2σ + 1) *
        ∫_ξ₁^ξ exp(real(c) * η^2) * η^(-d + 1 / σ) / (η^(1 - d) * exp(real(c) * η^2)) * η^((v - 1 / σ) * (2σ + 1)) dη
```
The last integral can be simplified to
```
∫_ξ₁^ξ η^(2σ * v + v - 3) dη = -(ξ₁^(2σ * v + v - 2) - ξ^(2σ * v + v - 2)) / (2σ * v + v - 2)
```
We can thus take
```
C₁ = -abs(1 + im * δ) / abs(1 - im * ϵ) * CE / abs(2c) *
    exp(imag(b - a) * π) * abs(c^b) * norm(u, v)^(2σ + 1) /
    (2σ * v + v - 2)
```

# Bounding `abs(g(ξ))`
Focusing on the integral we have
```
abs(∫_ξ^∞ P(η) / W(η) * abs(u(η))^2σ * u(η) dη) <=
    ∫_ξ^∞ abs(P(η)) / abs(W(η)) * abs(u(η))^(2σ + 1) dη
```
Using that
```
abs(P(η)) <= CP * η^(-1 / σ)
abs(W(η)) =  abs(2c) * exp(-imag(b - a) * π) * abs(c^-b) * η^(1 - d) * exp(real(c) * η^2)
abs(u(η)) <= norm(u, v) * η^(v - 1 / σ)
```
we get
```
abs(∫_ξ^∞ P(η) / W(η) * abs(u(η))^2σ * u(η) dη) <=
    CP / abs(2c) * exp(imag(b - a) * π) * abs(c^b) * norm(u, v)^(2σ + 1) *
        ∫_ξ^∞ η^(-1 / σ) / (η^(1 - d) * exp(real(c) * η^2)) * η^((v - 1 / σ) * (2σ + 1)) dη
```
The last integral can be simplified to
```
∫_ξ^∞ η^(-2 / σ + 2σ * v + v + d - 3) * exp(-real(c) * η^2) dη
```
By factoring out `exp(-real(c) * ξ^2)` and using that `real(c) = κ * ϵ
/ (1 + ϵ)^2 >= 0` we get
```
∫_ξ^∞ η^(-2 / σ + 2σ * v + v + d - 3) * exp(-real(c) * η^2) dη =
    exp(-real(c) * ξ^2) * ∫_ξ^∞ η^(-2 / σ + 2σ * v + v + d - 3) * exp(-real(c) * (η^2 - ξ^2)) dη <=
    exp(-real(c) * ξ^2) * ∫_ξ^∞ η^(-2 / σ + 2σ * v + v + d - 3) dη =
    exp(-real(c) * ξ^2) * η^(-2 / σ + 2σ * v + v d - 2) / (2 / σ - 2σ * v - v - d + 3)
```
We can thus take
```
C₂ = abs(1 + im * δ) / abs(1 - im * ϵ) * CP / abs(2c) *
    exp(imag(b - a) * π) * abs(c^b) * norm(u, v)^(2σ + 1) /
    (2 / σ - 2σ * v - v - d + 3)
```

# Bounding `norm(u, v)`
The bound for `norm(u, v)` is based on finding a ball that contains
the fix point of `T`.  If we let
```
CT₁, CT₂ = GinzburgLandauSelfSimilarSingular.C_T1(v, κ, p, ξ₁)
CT₃ = GinzburgLandauSelfSimilarSingular.C_T2(v, κ, p, ξ₁)
```
then if
```
CT₁ * r₁ * ξ₁^-v + CT₂ * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) <= ρ
```
and
```
2CT₃ * ρ^2σ * ξ₁^(-2 + 2σ * v) < 1
```
there is a fix point in the ball with radius `ρ`, hence `norm(u, v) <
ρ` in that case. We want to find the minimum `ρ` such that both above
inequalities are satisfied.

The second inequality is satisfied whenever
```
ρ < (ξ₁^(2 - 2σ * v) / 2CT₃)^(1 / 2σ)
```

For the first inequality we need
```
ρ * (1 - CT₂ * ξ₁^(-2 + 2σ * v) * ρ^2σ) >= CT₁ * r₁ * ξ₁^-v
```
The left hand side is clearly positive for
```
0 < ρ < ξ₁^(2 - 2σ * v) / CT₂
```
We want to pick `ρ` to be the minimum value so that the left hand side
is greater than `CT₁ * r₁ * ξ₁^-v`. We achieve this by solving the
equation
```
ρ * (1 - CT₂ * ξ₁^(-2 + 2σ * v) * ρ^2σ) = CT₁ * r₁ * ξ₁^-v
```
"""
function C_fix_point(r₁::Arb, v::Arb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)
    d, ω, σ, ϵ, δ = p.d, p.ω, p.σ, p.ϵ, p.δ

    a = Acb(1 / σ, ω / κ) / 2
    b = Arb(d // 2)
    c = Acb(0, -κ) / 2Acb(1, -ϵ)

    CP = C_P(κ, p, ξ₁)
    CE = C_E(κ, p, ξ₁)

    CT₁, CT₂ = GinzburgLandauSelfSimilarSingular.C_T1(v, κ, p, ξ₁)
    CT₃ = GinzburgLandauSelfSimilarSingular.C_T2(v, κ, p, ξ₁)

    ρ_bound = (ξ₁^(2 - 2σ * v) / 2CT₃)^(1 / 2σ)

    lower = Arf(1e-5)
    upper = lbound(min(ρ_bound, ξ₁^(2 - 2σ * v) / CT₂))

    roots, flags = ArbExtras.isolate_roots(lower, upper) do ρ
        ρ * (1 - CT₂ * ξ₁^(-2 + 2σ * v) * ρ^2σ) - CT₁ * r₁ * ξ₁^-v
    end

    @assert only(flags)

    root = ArbExtras.refine_root(Arb(only(roots))) do ρ
        ρ * (1 - CT₂ * ξ₁^(-2 + 2σ * v) * ρ^2σ) - CT₁ * r₁ * ξ₁^-v
    end

    # Take ρ slightly larger than the root so that we can verify the inequalities
    ρ = 1.01 * root

    @assert CT₁ * r₁ * ξ₁^-v + CT₂ * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) <= ρ
    @assert 2CT₃ * ρ^2σ * ξ₁^(-2 + 2σ * v) < 1

    C₁ =
        -abs(Acb(1, δ)) / abs(Acb(1, -ϵ)) * CE / abs(2c) *
        exp(imag(b - a) * π) *
        abs(c^b) *
        ρ^(2σ + 1) / (2σ * v + v - 2)

    C₂ =
        abs(Acb(1, δ)) / abs(Acb(1, -ϵ)) * CP / abs(2c) *
        exp(imag(b - a) * π) *
        abs(c^b) *
        ρ^(2σ + 1) / (2 / σ - 2σ * v - v - d + 3)

    return C₁, C₂
end

"""
    C_fix_point_dξ(r₁::Arb, v::Arb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)

Return `C₃` and `C₄` such that the fix point `u(ξ)` of the operator
`T` satisfies
```
u'(ξ) = γ * P'(ξ) + P'(ξ) * f(ξ) + P(ξ) * f'(ξ) + E'(ξ) * g(ξ) + E(ξ) * g'(ξ)
```
with
```
abs(f'(ξ)) <= C₃ * ξ^(2σ * v + v - 3)
abs(g'(ξ)) <= C₄ * exp(-real(c) * ξ^2) * ξ^(-2 / σ + 2σ * v + v + d - 3)
```
and `abs(f(ξ))` and `abs(g(ξ))` bounded as in [`C_fix_point`](@ref).
This holds for for `ξ >= ξ₁` and `abs(γ) <= r₁`. In particular for `ξ
= ξ₁` we have
```
u(ξ) = γ * P(ξ) + P(ξ) * f'(ξ) + E'(ξ) * g(ξ) + E(ξ) * g'(ξ)
```
Here
```
real(c) = κ * ϵ / (1 + ϵ)^2
```

# Formulas for `f'` and `g'`
As in [`C_fix_point`](@ref) we have
```
f(ξ) = -(1 + im * δ) / (1 - im * ϵ) * ∫_ξ₁^ξ E(η) / W(η) * abs(u(η))^2σ * u(η) dη

g(ξ) = -(1 + im * δ) / (1 - im * ϵ) * ∫_ξ^∞ P(η) / W(η) * abs(u(η))^2σ * u(η) dη
```
Differentiating we directly get
```
f'(ξ) = -(1 + im * δ) / (1 - im * ϵ) * E(ξ) / W(ξ) * abs(u(ξ))^2σ * u(ξ) dη

g'(ξ) = (1 + im * δ) / (1 - im * ϵ) * P(ξ) / W(ξ) * abs(u(ξ))^2σ * u(ξ) dη
```

# Bounding `abs(f'(ξ))`
Using that
```
abs(E(ξ)) <= CE * exp(real(c) * ξ^2) * ξ^(-d + 1 / σ)
abs(W(ξ)) =  abs(2c) * exp(-imag(b - a) * π) * abs(c^-b) * ξ^(1 - d) * exp(real(c) * ξ^2)
abs(u(ξ)) <= norm(u, v) * ξ^(v - 1 / σ)
```
we get
```
abs(f'(ξ)) <= abs(1 + im * δ) / abs(1 - im * ϵ) *
    CE / abs(2c) * exp(imag(b - a) * π) * abs(c^b) * norm(u, v)^(2σ + 1) *
    ξ^(-d + 1 / σ) / ξ^(1 - d) * ξ^((v - 1 / σ) * (2σ + 1))
```
We have
```
ξ^(-d + 1 / σ) / ξ^(1 - d) * ξ^((v - 1 / σ) * (2σ + 1)) = ξ^(2σ * v + v - 3)
```
and can take
```
C₃ = abs(1 + im * δ) / abs(1 - im * ϵ) * CE / abs(2c) *
    exp(imag(b - a) * π) * abs(c^b) * norm(u, v)^(2σ + 1)
```

# Bounding `abs(g(ξ))`
Using that
```
abs(P(ξ)) <= CP * ξ^(-1 / σ)
abs(W(ξ)) =  abs(2c) * exp(-imag(b - a) * π) * abs(c^-b) * ξ^(1 - d) * exp(real(c) * ξ^2)
abs(u(ξ)) <= norm(u, v) * ξ^(v - 1 / σ)
```
we get
```
abs(g'(ξ)) <= abs(1 + im * δ) / abs(1 - im * ϵ) *
    CP / abs(2c) * exp(imag(b - a) * π) * abs(c^b) * norm(u, v)^(2σ + 1) *
    ξ^(-1 / σ) / ξ^(1 - d) * exp(-real(c) * ξ^2) * ξ^((v - 1 / σ) * (2σ + 1))
```
We have
```
ξ^(-1 / σ) / ξ^(1 - d) * exp(-real(c) * ξ^2) * ξ^((v - 1 / σ) * (2σ + 1)) =
    exp(-real(c)) * ξ^(-2 / σ + 2σ * v + v + d - 3)
```
and can take
```
C₄ = abs(1 + im * δ) / abs(1 - im * ϵ) * CP / abs(2c) *
    exp(imag(b - a) * π) * abs(c^b) * norm(u, v)^(2σ + 1)
```

# Bounding `norm(u, v)`
This is done as in [`C_fix_point`](@ref).
"""
function C_fix_point_dξ(r₁::Arb, v::Arb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)
    d, ω, σ, ϵ, δ = p.d, p.ω, p.σ, p.ϵ, p.δ

    a = Acb(1 / σ, ω / κ) / 2
    b = Arb(d // 2)
    c = Acb(0, -κ) / 2Acb(1, -ϵ)

    CP = C_P(κ, p, ξ₁)
    CE = C_E(κ, p, ξ₁)

    CT₁, CT₂ = GinzburgLandauSelfSimilarSingular.C_T1(v, κ, p, ξ₁)
    CT₃ = GinzburgLandauSelfSimilarSingular.C_T2(v, κ, p, ξ₁)

    ρ_bound = (ξ₁^(2 - 2σ * v) / 2CT₃)^(1 / 2σ)

    lower = Arf(1e-5)
    upper = lbound(min(ρ_bound, ξ₁^(2 - 2σ * v) / CT₂))

    roots, flags = ArbExtras.isolate_roots(lower, upper) do ρ
        ρ * (1 - CT₂ * ξ₁^(-2 + 2σ * v) * ρ^2σ) - CT₁ * r₁ * ξ₁^-v
    end

    @assert only(flags)

    root = ArbExtras.refine_root(Arb(only(roots))) do ρ
        ρ * (1 - CT₂ * ξ₁^(-2 + 2σ * v) * ρ^2σ) - CT₁ * r₁ * ξ₁^-v
    end

    # Take ρ slightly larger than the root so that we can verify the
    # inequalities
    ρ = 1.01 * root

    @assert CT₁ * r₁ * ξ₁^-v + CT₂ * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) <= ρ
    @assert 2CT₃ * ρ^2σ * ξ₁^(-2 + 2σ * v) < 1

    C₃ =
        abs(Acb(1, δ)) / abs(Acb(1, -ϵ)) * CE / abs(2c) *
        exp(imag(b - a) * π) *
        abs(c^b) *
        ρ^(2σ + 1)

    C₄ =
        abs(Acb(1, δ)) / abs(Acb(1, -ϵ)) * CP / abs(2c) *
        exp(imag(b - a) * π) *
        abs(c^b) *
        ρ^(2σ + 1)

    return C₃, C₄
end

"""
    C_fix_point_dγ(r₁::Arb, v::Arb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)

Return `C` such that the fix point `u(ξ)` of the operator
`T` satisfies
```
diff(u(ξ), γ) = P(ξ) + h(ξ)
```
with
```
abs(h(ξ)) <= C * ξ^XXX
```
Here we use `diff(u(ξ₁), γ)` to denote the derivative of `u(ξ₁)` with
respect to `γ`. This holds for `ξ >= ξ₁` and `abs(γ) <= r₁`.

# Splitting `T`
Since `u` is a fix point of `T` we have
```
u(ξ) = γ * P(ξ) - (1 + im * δ) * ∫_ξ₁^∞ K(ξ, η) * abs(u(η))^2σ * u(η) dη
```
Differentiating with respect to `γ` we get
```
diff(u(ξ), γ) = P(ξ) - (1 + im * δ) * ∫_ξ₁^∞ K(ξ, η) * diff(abs(u(η))^2σ * u(η), γ) dη
```
We let
```
h(ξ) = (1 + im * δ) * ∫_ξ₁^∞ K(ξ, η) * diff(abs(u(η))^2σ * u(η), γ) dη
```

# Simplifying `diff(abs(u(η))^2σ * u(η), γ)`
We have
```
diff(abs(u(η))^2σ * u(η), γ) = diff(abs(u(η))^2σ, γ) * u(η) + abs(u(η))^2σ * diff(u(η), γ)
```
Furthermore if we let `u(η) = x(η) + im * y(η)` we have
```
diff(abs(u(η))^2σ, γ) = diff((x(η)^2 + y(η)^2)^σ, γ)
    = σ * (x(η)^2 + y(η)^2)^(σ - 1) * diff(x(η)^2 + y(η)^2, γ)
    = 2σ * abs(u(η))^(2σ - 2) * (x(η) * diff(x(η), γ) + y(η) * diff(y, γ))
```
If note that
```
conj(u(η)) * diff(u(η), γ) + u(η) * conj(diff(u(η), γ)) = x(η) * diff(x(η), γ) + y(η) * diff(y, γ)
```
Hence
```
diff(abs(u(η))^2σ, γ) = 2σ * abs(u(η))^(2σ - 2) * (
    conj(u(η)) * diff(u(η), γ) + u(η) * conj(diff(u(η), γ))
)
```
This gives us
```
diff(abs(u(η))^2σ * u(η), γ) =
    2σ * abs(u(η))^(2σ - 2) * conj(u(η)) * diff(u(η), γ) +
    2σ * abs(u(η))^(2σ - 2) * u(η) * conj(diff(u(η), γ)) +
    abs(u(η))^2σ * diff(u(η), γ)
```

# Bounding `abs(h(ξ))`
We have
```
abs(h(ξ)) <= abs(1 + im * δ) * ∫_ξ₁^∞ abs(K(ξ, η)) * abs(diff(abs(u(η))^2σ * u(η), γ)) dη
```
With the expression for `diff(abs(u(η))^2σ * u(η), γ)` from above we
get
```
abs(diff(abs(u(η))^2σ * u(η), γ)) <=
    2σ * abs(u(η))^(2σ - 1) * abs(diff(u(η), γ)) +
    2σ * abs(u(η))^(2σ - 1) * abs(diff(u(η), γ)) +
    abs(u(η))^2σ * abs(diff(u(η), γ))
= abs(u(η))^(2σ - 1) * (2σ + abs(u(η))) * abs(diff(u(η), γ))
```
This gives us
```
abs(h(ξ)) <= abs(1 + im * δ) * ∫_ξ₁^∞ abs(K(ξ, η)) *
    abs(u(η))^(2σ - 1) * (2σ + abs(u(η))) * abs(diff(u(η), γ)) dη
```
Using that
```
abs(u(η)) <= norm(u, v) * η^(v - 1 / σ)
abs(diff(u(η), γ)) <= norm(diff(u, γ), v) * η^(v - 1 / σ)
```
we get
```
abs(h(ξ)) <= abs(1 + im * δ) * norm(u, v)^(2σ - 1) * norm(diff(u, γ), v) *
    ∫_ξ₁^∞ abs(K(ξ, η)) * (2σ + norm(u, v) * η^(v - 1 / σ)) * η^((v - 1 / σ) * 2σ) dη
```

"""
function C_fix_point_dγ(r₁::Arb, v::Arb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)
    d, ω, σ, ϵ, δ = p.d, p.ω, p.σ, p.ϵ, p.δ
    q = 2
    a = Acb(1 / σ, ω / κ) / 2
    b = Arb(d // 2)
    c = Acb(0, -κ) / 2Acb(1, -ϵ)

    CP = C_P(κ, p, ξ₁)
    CE = C_E(κ, p, ξ₁)

    # TODO
    C = one(CP)

    return C
end

"""
    solution_infinity_γ(κ::Arb, μ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb, v::Arb; u_0::Acb)

Compute enclosure of `γ` occurring in the fixed point equation. If the
fixed point cannot be proved to exist then returns an indeterminate
value.

**IMPROVE:** Add more documentation
"""
function solution_infinity_γ(
    κ::Arb,
    μ::Arb,
    p::AbstractGLParams{Arb},
    ξ₁::Arb,
    v::Arb,
    u_0::Acb,
)
    # Compute approximate γ
    γ₀ = u_0 / P(ξ₁, (p, κ))

    # Enclose u_T(ξ₁, (p, κ, γ₀))
    u_T = let
        _, C₂ = C_fix_point(abs(γ₀), v, κ, p, ξ₁)

        real_c = κ * p.ϵ / (1 + p.ϵ)^2

        g_bound = add_error(
            zero(C₂),
            C₂ * exp(-real_c * ξ₁^2) * ξ₁^(-2 / p.σ + 2p.σ * v + v + p.d - 2),
        )

        γ₀ * P(ξ₁, (p, κ)) + E(ξ₁, (p, κ)) * g_bound
    end

    # TODO: Choose the radius dynamically
    γ_ball = add_error(γ₀, Mag(1e-3))

    u_T_dγ = let
        C₃ = C_fix_point_dγ(abs(γ_ball), v, κ, p, ξ₁)

        # FIXME: This is not the correct expression
        real_c = κ * p.ϵ / (1 + p.ϵ)^2

        h_bound = add_error(
            zero(C₃),
            C₃ *
            abs(E(ξ₁, (p, κ))) *
            exp(-real_c * ξ₁^2) *
            ξ₁^(-2 / p.σ + 2p.σ * v + v + p.d - 2),
        )

        P(ξ₁, (p, κ)) + h_bound
    end

    # Perform one Newton iteration
    γ_ball_new = γ₀ - (u_T - u_0) / u_T_dγ

    success = Arblib.contains_interior(γ_ball, γ_ball_new)

    if success
        return γ_ball_new
    else
        return indeterminate(γ_ball_new)
    end
end
