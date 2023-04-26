"""
    C_hypgeom_u(a::Acb, b::Acb, z₁::Arb)

Return `C` such that
```
abs(hypgeom_u(a, b, z)) <= C * abs(z^-a)
```
for `z` such that `abs(imag(z)) > abs(imag(z₁))` and `abs(z) > abs(z₁)`

It makes use of the asymptotic expansion given by
https://fungrim.org/entry/d1b3b5/ together with error bounds for the
remainder term from https://fungrim.org/entry/461a54/.

The asymptotic expansion with `n` terms is given by
```
z^-a * sum(0:n-1) do k
    rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-z)^k)
end
```
With the absolute value of the remainder bounded by
```
R / abs(z)^(a + n)
```
where
```
R = abs(rising(a, n) * rising(a - b + 1, n) / factorial(n)) *
        2 / (1 - s) * exp(π * ρ / ((1 - s) * abs(z₁)))
```
with
```
s = abs(b - 2a) / abs(z₁)
ρ = abs(a^2 - a * b + b / 2) + s * (1 + s / 4) / (1 - s)^2
```

We get a bound for the asymptotic expansion, skipping the `z^-a`
factor, as
```
abs(
    sum(0:n-1) do k
        rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-z)^k)
    end
) <=
sum(0:n-1) do k
    abs(rising(a, k) * rising(a - b + 1, k)) / (factorial(k) * abs(z)^k)
end <=
sum(0:n-1) do k
    abs(rising(a, k) * rising(a - b + 1, k)) / (factorial(k) * abs(z₁)^k)
end
```
For the remainder term we get, also skipping the `z^-a` factor,
```
R / abs(z)^n <= R / abs(z₁)^n
```
"""
function C_hypgeom_u(a::Acb, b::Acb, z₁::Acb, n::Integer = 5)
    abs(imag(z₁)) > abs(imag(b - 2a)) ||
        throw(ArgumentError("assuming abs(imag(z)) > abs(imag(z₁))"))

    C = sum(0:n-1) do k
        abs(rising(a, k) * rising(a - b + 1, k)) / (factorial(k) * abs(z₁)^k)
    end

    s = abs(b - 2a) / abs(z₁)
    ρ = abs(a^2 - a * b + b / 2) + s * (1 + s / 4) / (1 - s)^2

    R =
        abs(rising(a, n) * rising(a - b + 1, n) / factorial(n)) * 2 / (1 - s) *
        exp(π * ρ / ((1 - s) * abs(z₁)))

    C += R / abs(z₁)^n

    return C
end

"""
    C_P(κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)

Return `C` such that
```
abs(P(ξ, (p, κ))) <= C * ξ^(-1 / σ)
```
for `ξ >= ξ₁ > 0`.

Let
```
a = (1 / σ + im * ω / κ) / 2
b = d / 2
c = -im * κ / 2(1 - im * ϵ)
```
We then have
```
P(ξ) = hypgeom_u(a, b, c * ξ^2)
```

If we let `CU = C_hypgeom_u(a, b, c * ξ₁)` then we have for `ξ` such
that `abs(imag(c * ξ)) > abs(imag(c * ξ₁))` and `abs(c * ξ₁) > abs(c *
ξ₁)` that
```
abs(hypgeom_u(a, b, c * ξ)) <= CU * abs((c * ξ)^-a)
    = CU * abs(c^-a) * abs(ξ^-a)
    = CU * abs(c^-a) * ξ^(-1 / σ)
```
In particular we have that if `ξ >= ξ₁` then both the requirements are
satisfied. We can thus take
```
C = CU * abs(c^-a)
```
"""
function C_P(κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)
    ξ₁ > 0 || throw(ArgumentError("assuming ξ₁ > 0"))

    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = Acb(1 / σ, ω / κ) / 2
    b = Acb(d // 2)
    c = Acb(0, -κ) / 2Acb(1, -ϵ)

    C = C_hypgeom_u(a, b, c * ξ₁^2) * abs(c^-a)

    return C
end

"""
    C_E(κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)

Return `C` such that
```
abs(E(ξ, (p, κ))) <= C * exp(real(c) * ξ^2) * ξ^(-d + 1 / σ)
```
for `ξ >= ξ₁ > 0`, with `c` as below. Note that if `ϵ = 0` then
`real(z) = 0` and the exponential term can be removed.

Se also [`C_P`](@ref).

Let
```
a = (1 / σ + im * ω / κ) / 2
b = d / 2
c = -im * κ / 2(1 - im * ϵ)
```
We then have
```
E(ξ) = exp(c * ξ^2) * hypgeom_u(b - a, b, -c * ξ^2)
```

If we let `CU = C_hypgeom_u(b - a, b, -c * ξ₁)` then we have for `ξ`
such that `abs(imag(c * ξ)) > abs(imag(c * ξ₁))` and `abs(c * ξ₁) >
abs(c * ξ₁)` that
```
abs(hypgeom_u(b - a, b, -c * ξ)) <= CU * abs((-c * ξ^2)^(a - b))
    = CU * abs((-c)^(a - b)) * abs(ξ^(2(a - b)))
    = CU * abs((-c)^(a - b)) * ξ^(-d - 1 / σ)
```
In particular we have that if `ξ >= ξ₁` then `z` satisfies both the
requirements. If we let `C = CU * abs((-c)^(a - b))` and use that
`abs(exp(c * ξ^2)) = exp(real(c) * ξ^2)` we thus get
```
E(ξ) <= C * exp(real(c) * ξ^2) * ξ^(-d - 1 / σ)
```
as required.
"""
function C_E(κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)
    ξ₁ > 0 || throw(ArgumentError("assuming ξ₁ > 0"))

    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = Acb(1 / σ, ω / κ) / 2
    b = Acb(d // 2)
    c = Acb(0, -κ) / 2Acb(1, -ϵ)

    C = C_hypgeom_u(b - a, b, c * ξ₁^2) * abs((-c)^(a - b))

    return C
end

"""
    C_K(κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)

Return `C` such that
```
abs(K(ξ, η, (p, κ))) <= C * ξ^(-1 / σ) * η^(1 / σ - 1)
```
for `ξ₁ <= η <= ξ` and
```
abs(K(ξ, η, (p, κ))) <= C * ξ^(-d + 1 / σ) * η^(-1 - 1 / σ + d)
```
for `ξ₁ <= ξ <= η`.

For we have `η <= ξ`
```
K(ξ, η) = -1 / (1 - im * ϵ) * P(ξ) * E(η) / W(η)
```
and for `ξ <= η`
```
K(ξ, η) = -1 / (1 - im * ϵ) * E(ξ) * P(η) / W(η)
```

# Handling `W(ξ)`
Let
```
a = (1 / σ + im * ω / κ) / 2
b = d / 2
c = -im * κ / 2(1 - im * ϵ)
```
We then have
```
W(ξ) = 2c * exp(im * (b - a) * π) * ξ * (c * ξ^2)^-b * exp(c * ξ^2)
```
For the absolute value we get
```
abs(W(ξ)) = abs(2c) * exp(-imag(b - a) * π) * ξ * abs(c * ξ^2)^-b * exp(real(c) * ξ^2)
```
From the definition of `b` we can write this as
```
abs(W(ξ)) = abs(2c) * exp(-imag(b - a) * π) * abs(c^-b) * ξ^(1 - d) * exp(real(c) * ξ^2)
```

# Bound for `η <= ξ`
Let `CP` and `CE` be the coefficients from [`C_P`](@ref) and
[`C_E`](@ref), we then have
```
abs(P(ξ)) <= CP * ξ^(-1 / σ)
abs(E(η)) <= CE * exp(real(c) * ξ^2) * η^(-d + 1 / σ)
```
This gives us
```
abs(K(ξ, η)) <= CP * CE / abs(1 - im * ϵ) * ξ^(-1 / σ) * η^(-d + 1 / σ) *
    exp(real(c) * η^2) /
    (
        abs(2c) * exp(-imag(b - a) * π) * abs(c^-b) * η^(1 - d) * exp(real(c) * η^2)
    )
= CP * CE * exp(imag(b - a) * π) * abs(c^b) / (abs(2c) * abs(1 - im * ϵ)) *
    ξ^(-1 / σ) * η^(-1 + 1 / σ)
```
Note that `abs(2c) * abs(1 - im * ϵ) = κ`, we can thus take
```
C = CP * CE * abs(c^-b) * exp(imag(b - a) * π) / κ
```

# Bound for `ξ <= η`
In this case we get
```
abs(K(ξ, η)) <= CP * CE / abs(1 - im * ϵ) * ξ^(-d + 1 / σ) * η^(-1 / σ) *
    exp(real(c) * ξ^2) /
    (
        abs(2c) * exp(-imag(b - a) * π) * abs(c^-b) * η^(1 - d) * exp(real(c) * η^2)
    )
= CP * CE * exp(imag(b - a) * π) * abs(c^b) / (abs(2c) * abs(1 - im * ϵ)) *
    exp(real(c) * (ξ^2 - η^2)) * ξ^(-d + 1 / σ) * η^(-1 - 1 / σ + d)
```
Compared to the previous case the only difference is
```
exp(real(c) * (ξ^2 - η^2))
```
Since
```
real(c) = κ * ϵ / (1 + ϵ)^2 >= 0
```
and `ξ^2 - η^2 < 0` this is bounded by `1` and the same `C` as above
gives an upper bound.
"""
function C_K(κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    κ >= 0 || throw(ArgumentError("assuming κ >= 0"))
    ϵ >= 0 || throw(ArgumentError("assuming ϵ >= 0"))

    a = Acb(1 / σ, ω / κ) / 2
    b = Arb(d // 2)
    c = Acb(0, -κ) / 2Acb(1, -ϵ)

    CP = C_P(κ, p, ξ₁)
    CE = C_E(κ, p, ξ₁)

    C = CP * CE * exp(imag(b - a) * π) * abs(c^b) / κ

    return C
end

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
- **TODO:** Check that this is enough
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
