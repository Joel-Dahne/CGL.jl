"""
    C_P(κ::Arb, p::AbstractGLParams{Arb})

Return `C` such that
```
abs(P(ξ, (p, κ))) <= C * ξ^(-1 / σ)
```
for `ξ >= 1`.

Let
```
a = (1 / σ + im * ω / κ) / 2
b = d / 2
z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
```
We then have
```
P(ξ) = hypgeom_u(a, b, z)
```
Assuming that `real(a) > 0` and `-π < angle(z) < π` we have the
integral representation
```
hypgeom_u(a, b, z) = z^(-a) / gamma(a) * ∫ exp(-s) * s^(a - 1) * (1 + s / z)^(b - a - 1) ds
```
with the integration taken from `0` to `∞`. Note that since `σ, ω, κ`
are all real we have `real(a) > 0` as long as `σ > 0`. Furthermore,
since `ξ >= 1` the requirement for `z` is also satisfied.

This gives us
```
abs(hypgeom_u(a, b, z)) <= abs(z)^(-real(a)) / abs(gamma(a)) *
    ∫ exp(-s) * s^(real(a) - 1) * (1 + abs(s / z))^real(b - a - 1) ds
```

We have
```
abs(z)^(-real(a)) = abs(-im * κ / (1 - im * ϵ) / 2)^(-real(a)) * ξ^(-2real(a))
```

Let
```
I = ∫ exp(-s) * s^(real(a) - 1) * (1 + abs(s / z))^real(b - a - 1) ds
```
Assuming that `real(b - a - 1) <= 0` we have `(1 + abs(s / z))^real(b -
a - 1) <= 1` and hence
```
I <= ∫ exp(-s) * s^(real(a) - 1) ds = gamma(real(a))
```
If `real(b - a - 1) <= 1` we have
```
(1 + abs(s / z))^real(b - a - 1) <= 1 + abs(s / z) * (1 + abs(s / z))^real(b - a - 2)
    <= 1 + abs(s / z)
```
and hence
```
I <= ∫ exp(-s) * s^(real(a) - 1) ds + 1 / abs(z) * ∫ exp(-s) * s^real(a) ds
  = gamma(real(a)) + gamma(real(a + 1)) / abs(z)
```
Where we can use that
```
1 / abs(z) = inv(κ / abs(1 - im * ϵ) * ξ^2 / 2)
           = 2abs(1 - im * ϵ) / (κ * ξ^2)
           <= 2abs(1 - im * ϵ) / κ
```
- **TODO:** Handle also `b - a - 1 > 1`.

Combining this we get
```
abs(P(ξ)) <= abs(-im * κ / (1 - im * ϵ) / 2)^(-real(a)) * I / abs(gamma(a)) * ξ^(-2real(a))
```
If we let
```
C = abs(-im * κ / (1 - im * ϵ) / 2)^(-real(a)) * I / abs(gamma(a))
```
and use that `real(a) = 1 / 2σ` we get
```
abs(P(ξ)) <= C * ξ^(-1 / σ)
```

**IMPROVE:** This gives fairly poor bounds and can be improved. In
particular if we restrict it to `ξ >= ξ₁` for some `ξ₁` we can get a
much better bound for large values.
"""
function C_P(κ::Arb, p::AbstractGLParams{Arb})
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    σ > 0 || throw(ArgumentError("assuming σ > 0"))

    a = Acb(1 / σ, ω / κ) / 2
    b = Arb(d // 2)

    if real(b - a - 1) <= 0
        I = gamma(real(a))
    elseif real(b - a - 1) <= 1
        I = gamma(real(a)) + gamma(real(a + 1)) * 2abs(1 - im * ϵ) / κ
    else
        error("currently only implements b - a - 1 <= 1")
    end

    C = abs(-im * κ / (1 - im * ϵ) / 2)^(-real(a)) * I / abs(gamma(a))

    return C
end

"""
    C_E(κ::Arb, p::AbstractGLParams{Arb})

Return `C` such that
```
abs(E(ξ, (p, κ))) <= C * exp(real(z)) * ξ^(-d + 1 / σ)
```
for `ξ >= 1`, with `z` as below. Note that if `ϵ = 0` then `real(z) =
0` and the exponential term can be removed.

Se also [`C_E`](@ref).

Let
```
a = (1 / σ + im * ω / κ) / 2
b = d / 2
z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
```
We then have
```
E(ξ) = exp(z) * hypgeom_u(b - a, b, -z)
```
Assuming that `real(b - a) > 0` and `-π < angle(-z) < π` we have the
integral representation
```
E(ξ) = -exp(z) * z^(-(b - a)) / gamma(b - a) * ∫ exp(-s) * s^(b - a - 1) * (1 - s / z)^(a - 1) ds
```
with the integration taken from `0` to `∞`. Since `ξ >= 1` the
requirement for `z` is always satisfied, the requirement for `real(b -
a)` we have to check

This gives us
```
abs(E(ξ)) <= exp(real(z)) * abs(z)^(-real(b - a)) / abs(gamma(b - a)) *
    ∫ exp(-s) * s^real(b - a - 1) * (1 + s / z)^real(a - 1) ds
```

We have
```
abs(z)^(-real(b - a)) = abs(-im * κ / (1 - im * ϵ) / 2)^(-real(b - a)) * ξ^(-2real(b - a))
```

Let
```
I = ∫ exp(-s) * s^real(b - a - 1) * (1 + s / z)^real(a - 1) ds
```
Assuming that `real(a - 1) <= 0` we have `(1 + abs(s / z))^real(a - 1) <= 1` and hence
```
I <= ∫ exp(-s) * s^real(b - a - 1) ds = gamma(real(b - a))
```
- **TODO:** Do we need to also handle `b - a - 1 > 0`.

Combining this we get
```
abs(E(ξ)) <= abs(-im * κ / (1 - im * ϵ) / 2)^(-real(b - a)) * I / abs(gamma(b - a)) * ξ^(-2real(b - a))
```
If we let
```
C = abs(-im * κ / (1 - im * ϵ) / 2)^(-real(b - a)) * I / abs(gamma(b - a))
```
and use that `real(b - a) = d / 2 - 1 / 2σ` we get
```
abs(E(ξ)) <= C * ξ^(-d - 1 / σ)
```

**IMPROVE:** This gives fairly poor bounds and can be improved. In
particular if we restrict it to `ξ >= ξ₁` for some `ξ₁` we can get a
much better bound for large values.
"""
function C_E(κ::Arb, p::AbstractGLParams{Arb})
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = Acb(1 / σ, ω / κ) / 2
    b = Arb(d // 2)

    real(b - a) > 0 || error("assuming b - a > 0")

    if real(a - 1) <= 0
        I = gamma(real(b - a))
    else
        error("currently only implements a - 1 <= 0")
    end

    C = abs(-im * κ / (1 - im * ϵ) / 2)^(-real(b - a)) * I / abs(gamma(b - a))

    return C
end

"""
    C_K(κ::Arb, p::AbstractGLParams{Arb})

Return `C` such that
```
abs(K(ξ, η, (p, κ))) <= C * ξ^(-1 / σ) * η^(1 / σ - 1)
```
for `1 <= η <= ξ` and
```
abs(K(ξ, η, (p, κ))) <= C * ξ^(-d + 1 / σ) * η^(-1 - 1 / σ + d)
```
for `1 <= ξ <= η`.

Let
```
a = (1 / σ + im * ω / κ) / 2
b = d / 2
z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
```
We then have for `1 <= η <= ξ`
```
K(ξ, η) = -1 / (1 - im * ϵ) * P(ξ) * E(η) / W(η)
```
and for `1 <= ξ <= η`
```
K(ξ, η) = -1 / (1 - im * ϵ) * E(ξ) * P(η) / W(η)
```

# Handling `W(ξ)`
In `K(ξ, η)` there is a division by
```
W(ξ) = -im * κ / (1 - im * ϵ) * exp(im * (b - a) * π) * ξ * z^-b * exp(z)
```
For which we have
```
abs(W(ξ)) = κ / abs(1 - im * ϵ) * exp(-imag(b - a) * π) * ξ *
    abs(z)^-b * exp(real(z))
```
From the definition of `z` and `b` we can write this as
```
abs(W(ξ)) = κ / abs(1 - im * ϵ) * exp(-imag(b - a) * π) * ξ^(1 - d) *
    abs(κ / 2(1 - im * ϵ))^-b * exp(real(z))
```

# Bound for `η <= ξ`
Let `CP` and `CE` be the coefficients from [`C_P`](@ref) and
[`C_E`](@ref), we then have
```
abs(P(ξ)) <= CP * ξ^(-1 / σ)
abs(E(η)) <= CE * exp(real(-im * κ / (1 - im * ϵ) * η^2 / 2)) * η^(-d + 1 / σ)
```
This gives us
```
abs(K(ξ, η)) <= CP * CE / abs(1 - im * ϵ) * ξ^(-1 / σ) * η^(-d + 1 / σ) *
    exp(real(-im * κ / (1 - im * ϵ) * η^2 / 2)) /
    (
        κ / abs(1 - im * ϵ) * exp(- π * imag(b - a)) * η^(1 - d) *
        abs(κ / (1 - im * ϵ))^-b * exp(real(-im * κ / (1 - im * ϵ) * η^2 / 2))
    )
= CP * CE * abs(κ / (1 - im * ϵ))^b * exp(- π * imag(b - a)) / κ *
    ξ^(-1 / σ) * η^(-1 + 1 / σ)
```
So we can take `C` to
```
C = CP * CE * abs(κ / (1 - im * ϵ))^b * exp(- π * imag(b - a)) / κ
```

# Bound for `ξ <= η`
In this case we get
```
abs(K(ξ, η)) <= CP * CE / abs(1 - im * ϵ) * ξ^(-d + 1 / σ) * η^(-1 / σ) *
    exp(real(-im * κ / (1 - im * ϵ) * ξ^2 / 2)) /
    (
        κ / abs(1 - im * ϵ) * exp(- π * imag(b - a)) * η^(1 - d) *
        abs(κ / (1 - im * ϵ))^-b * exp(real(-im * κ / (1 - im * ϵ) * η^2 / 2))
    )
= CP * CE * abs(κ / (1 - im * ϵ))^b * exp(- π * imag(b - a)) / κ *
    exp(real(-im * κ / (1 - im * ϵ) * ξ^2 / 2) - real(-im * κ / (1 - im * ϵ) * η^2 / 2)) *
    ξ^(-d + 1 / σ) * η^(-1 - 1 / σ + d)
```
Compared to the previous case the only difference is
```
exp(real(-im * κ / (1 - im * ϵ) * ξ^2 / 2) - real(-im * κ / (1 - im * ϵ) * η^2 / 2)) =
    exp(real(-im * κ / (1 - im * ϵ)) * (ξ^2 - η^2) / 2)
```
Since
```
real(-im * κ / (1 - im * ϵ)) = κ * ϵ / (1 + ϵ)^2 >= 0
```
and `ξ^2 - η^2 < 0` this is bounded by one and the same `C` gives an
upper bound.
"""
function C_K(κ::Arb, p::AbstractGLParams{Arb})
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    κ >= 0 || throw(ArgumentError("assuming κ >= 0"))
    ϵ >= 0 || throw(ArgumentError("assuming ϵ >= 0"))

    a = Acb(1 / σ, ω / κ) / 2
    b = Arb(d // 2)

    CP = C_P(κ, p)
    CE = C_E(κ, p)

    C = CP * CE * abs(κ / 2(1 - im * ϵ))^b * exp(imag(b - a) * π) / κ

    return C
end

"""
    C_T1(v::Arb, κ::Arb, p::AbstractGLParams{Arb})

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
function C_T1(v::Arb, κ::Arb, p::AbstractGLParams{Arb})
    d, ω, σ, ϵ, δ = p.d, p.ω, p.σ, p.ϵ, p.δ

    (2σ + 1)v < 2 + 2 / σ - d || throw(ArgumentError("assuming (2σ + 1)v < 2 + 2 / σ - d"))


    # TODO: This might follow from the other requirements?
    @assert -2 + (2σ + 1) * v < 0
    @assert -2 + (2σ + 1) * v + d - 2 / σ < 0

    C₁ = C_P(κ, p)

    C₂ =
        abs(1 + im * δ) *
        C_K(κ, p) *
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

TODO
"""
function C_T2(v::Arb, γ::Acb, κ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)
    return indeterminate(κ)
end
