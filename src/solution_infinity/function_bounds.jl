"""
    C_hypgeom_u(a::Acb, b::Acb, z₁::Acb)

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
R(z) / abs(z)^a
```
where
```
R(z) = abs(rising(a, n) * rising(a - b + 1, n) / (factorial(n) * z^n)) *
        2sqrt(1 + π * n / 2) / (1 - s) * exp(π * ρ / ((1 - s) * abs(z₁)))
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
end =
sum(0:n-1) do k
    abs(rising(a, k) * rising(a - b + 1, k)) / (factorial(k) * abs(z₁)^k)
end
```
For the remainder term it can also be checked that
```
R(z) <= R(z₁)
```
"""
function C_hypgeom_u(a::Acb, b::Acb, z₁::Acb, n::Integer = 5)
    abs(imag(z₁)) > abs(imag(b - 2a)) ||
        throw(ArgumentError("assuming abs(imag(z₁)) > abs(imag(b - 2a))"))

    C = sum(0:n-1) do k
        abs(rising(a, k) * rising(a - b + 1, k)) / (factorial(k) * abs(z₁)^k)
    end

    s = abs(b - 2a) / abs(z₁)
    ρ = abs(a^2 - a * b + b / 2) + s * (1 + s / 4) / (1 - s)^2

    R =
        abs(rising(a, n) * rising(a - b + 1, n) / (factorial(n) * abs(z₁)^n)) *
        2sqrt(1 + Arb(π) * n / 2) / (1 - s) * exp(π * ρ / ((1 - s) * abs(z₁)))

    C += R

    return C
end

"""
    C_hypgeom_u_dz(a::Acb, b::Acb, z₁::Acb)

Return `C` such that
```
abs(hypgeom_u_dz(a, b, z)) <= C * abs(z^(-a - 1))
```
for `z` such that `abs(imag(z)) > abs(imag(z₁))` and `abs(z) > abs(z₁)`

Uses that
```
hypgeom_u_dz(a, b, z) = -hypgeom_u(a + 1, b + 1, z) * a
```
"""
C_hypgeom_u_dz(a::Acb, b::Acb, z₁::Acb, n::Integer = 5) =
    abs(a) * C_hypgeom_u(a + 1, b + 1, z₁, n)

"""
    C_P(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)

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
function C_P(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    ξ₁ > 0 || throw(ArgumentError("assuming ξ₁ > 0"))

    a, b, c = _abc(κ, λ)

    C = C_hypgeom_u(a, b, c * ξ₁^2) * abs(c^-a)

    return C
end

"""
    C_E(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)

Return `C` such that
```
abs(E(ξ, (p, κ))) <= C * exp(real(c) * ξ^2) * ξ^(1 / σ - d)
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

If we let `CU = C_hypgeom_u(b - a, b, -c * ξ₁^2)` then we have for `ξ`
such that `abs(imag(c * ξ)) > abs(imag(c * ξ₁))` and `abs(c * ξ₁) >
abs(c * ξ₁)` that
```
abs(hypgeom_u(b - a, b, -c * ξ^2)) <= CU * abs((-c * ξ^2)^(a - b))
    = CU * abs((-c)^(a - b)) * abs(ξ^(2(a - b)))
    = CU * abs((-c)^(a - b)) * ξ^(1 / σ - d)
```
In particular we have that if `ξ >= ξ₁` then `z` satisfies both the
requirements. If we let `C = CU * abs((-c)^(a - b))` and use that
`abs(exp(c * ξ^2)) = exp(real(c) * ξ^2)` we thus get
```
E(ξ) <= C * exp(real(c) * ξ^2) * ξ^(1 / σ - d)
```
as required.
"""
function C_E(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    ξ₁ > 0 || throw(ArgumentError("assuming ξ₁ > 0"))

    a, b, c = _abc(κ, λ)

    C = C_hypgeom_u(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b))

    return C
end

# TODO
function C_P_dξ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    f = ξ -> ξ^(-1 / λ.σ - 1)

    # FIXME: This is only an approximation. It seems to be good
    # though.
    return 1.01 * abs(P_dξ(ξ₁, (λ, κ))) / f(ξ₁)
end

# TODO
function C_P_dξ_dξ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    f = ξ -> ξ^(-1 / λ.σ - 2)

    # FIXME: This is only an approximation. It seems to be good
    # though.
    return 1.01 * abs(P_dξ_dξ(ξ₁, (λ, κ))) / f(ξ₁)
end

# TODO
function C_E_dξ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    _, _, c = _abc(κ, λ)

    f = ξ -> exp(real(c) * ξ^2) * ξ^(1 / λ.σ - λ.d + 1)

    # FIXME: This is only an approximation. It seems to be good
    # though.
    return 1.01 * abs(E_dξ(ξ₁, (λ, κ))) / f(ξ₁)
end

function C_E_dξ_dξ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    _, _, c = _abc(κ, λ)

    f = ξ -> exp(real(c) * ξ^2) * ξ^(1 / λ.σ - λ.d + 2)

    # FIXME: This is only an approximation. It seems to be good
    # though.
    return 1.05 * abs(E_dξ_dξ(ξ₁, (λ, κ))) / f(ξ₁)
end

# TODO
function C_P_dκ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    f = ξ -> log(ξ) * ξ^(-1 / λ.σ)

    # FIXME: This is only an approximation. We multiply with 2 since
    # numerically it seems that the quotient is increasing in ξ.
    return 2 * abs(P_dκ(ξ₁, (λ, κ))) / f(ξ₁)
end

# TODO
function C_E_dκ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    _, _, c = _abc(κ, λ)

    f = ξ -> exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 2)

    # FIXME: This is only an approximation. It seems to be good
    # though.
    return 2 * abs(E_dκ(ξ₁, (λ, κ))) / f(ξ₁)
end

# TODO
function C_P_dξ_dκ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    f = ξ -> log(ξ) * ξ^(-1 / λ.σ - 1)

    # FIXME: This is only an approximation. We multiply with 3 since
    # numerically it seems that the quotient is increasing quite a lot
    # in ξ.
    return 3 * abs(P_dξ_dκ(ξ₁, (λ, κ))) / f(ξ₁)
end

# TODO: Might not need this
function C_E_dξ_dκ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    return indeterminate(ξ₁)
end

"""
    C_K(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)

Return `C` such that
```
abs(K(ξ, η, (p, κ))) <= C * ξ^(-1 / σ) * η^(1 / σ - 1)
```
for `ξ₁ <= η <= ξ` and
```
abs(K(ξ, η, (p, κ))) <= C * ξ^(1 / σ - d) * η^(d - 1 / σ - 1)
```
for `ξ₁ <= ξ <= η`.

Where for `η <= ξ` we have
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
abs(E(η)) <= CE * exp(real(c) * ξ^2) * η^(1 / σ - d)
```
This gives us
```
abs(K(ξ, η)) <= CP * CE / abs(1 - im * ϵ) * ξ^(-1 / σ) * η^(1 / σ - d) *
    exp(real(c) * η^2) /
    (
        abs(2c) * exp(-imag(b - a) * π) * abs(c^-b) * η^(1 - d) * exp(real(c) * η^2)
    )
= CP * CE * exp(imag(b - a) * π) * abs(c^b) / (abs(2c) * abs(1 - im * ϵ)) *
    ξ^(-1 / σ) * η^(1 / σ - 1)
```
Note that `abs(2c) * abs(1 - im * ϵ) = κ`, we can thus take
```
C = CP * CE * abs(c^-b) * exp(imag(b - a) * π) / κ
```

# Bound for `ξ <= η`
In this case we get
```
abs(K(ξ, η)) <= CP * CE / abs(1 - im * ϵ) * ξ^(1 / σ - d) * η^(-1 / σ) *
    exp(real(c) * ξ^2) /
    (
        abs(2c) * exp(-imag(b - a) * π) * abs(c^-b) * η^(1 - d) * exp(real(c) * η^2)
    )
= CP * CE * exp(imag(b - a) * π) * abs(c^b) / (abs(2c) * abs(1 - im * ϵ)) *
    exp(real(c) * (ξ^2 - η^2)) * ξ^(1 / σ - d) * η^(d - 1 / σ - 1)
```
Compared to the previous case the only difference is
```
exp(real(c) * (ξ^2 - η^2))
```
Since
```
real(c) = κ * ϵ / (1 + ϵ^2) >= 0
```
and `ξ^2 - η^2 < 0` this is bounded by `1` and the same `C` as above
gives an upper bound.
"""
function C_K(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    κ >= 0 || throw(ArgumentError("assuming κ >= 0"))
    λ.ϵ >= 0 || throw(ArgumentError("assuming ϵ >= 0"))

    a, b, c = _abc(κ, λ)

    CP = C_P(κ, λ, ξ₁)
    CE = C_E(κ, λ, ξ₁)

    sgn = sign(imag(c))

    C = CP * CE * exp(sgn * imag(b - a) * π) * abs(c^b) / κ

    return C
end

"""
    C_J_P(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})

Return `C` such that for
```
J_P(ξ) = -(1 + im * δ) / (1 - im * ϵ) * P(ξ) * inv(W(ξ))
```
we have
```
abs(J_P(ξ)) <= C * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d - 1)
```
for `ξ >= ξ₁`.
"""
function C_J_P(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    _, _, c = _abc(κ, λ)

    f = ξ -> exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d - 1)

    # FIXME: This is only an approximation
    return 1.01 * abs(J_P(ξ₁, (λ, κ))) / f(ξ₁)
end

"""
    C_J_E(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})

Return `C` such that for
```
J_E(ξ) = -(1 + im * δ) / (1 - im * ϵ) * E(ξ) * inv(W(ξ))
```
we have
```
abs(J_E(ξ)) <= C * ξ^(1 / σ - 1)
```
for `ξ >= ξ₁`.
"""
function C_J_E(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    f = ξ -> ξ^(1 / λ.σ - 1)

    # FIXME: This is only an approximation
    return 1.01 * abs(J_E(ξ₁, (λ, κ))) / f(ξ₁)
end

# TODO
function C_J_P_dξ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    _, _, c = _abc(κ, λ)

    f = ξ -> exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d)

    # FIXME: This is only an approximation. It seems to be good
    # though.
    return 1.01 * abs(J_P_dξ(ξ₁, (λ, κ))) / f(ξ₁)
end

# TODO
function C_J_E_dξ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    f = ξ -> ξ^(1 / λ.σ - 2)

    # FIXME: This is only an approximation. It seems to be good
    # though.
    return 1.01 * abs(J_E_dξ(ξ₁, (λ, κ))) / f(ξ₁)
end

"""
    C_J_P_dκ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})

Return `C` such that for `J_P_dκ(ξ)` being the derivative w.r.t. `κ`
of
```
J_P(ξ) = -(1 + im * δ) / (1 - im * ϵ) * P(ξ) * inv(W(ξ))
J_P_dκ(ξ) = -(1 + im * δ) / (1 - im * ϵ) * P_dκ(ξ) * inv(W(ξ))
```
we have
```
abs(J_P_dκ(ξ)) <= C * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d + 1)
```
for `ξ >= ξ₁`.
"""
function C_J_P_dκ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    _, _, c = _abc(κ, λ)

    f = ξ -> exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d + 1)

    # FIXME: This is only an approximation. It seems to be good
    # though.
    return 1.01 * abs(J_P(ξ₁, (λ, ArbSeries((κ, 1))))[1]) / f(ξ₁)
end

"""
    C_J_E_dκ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})

Return `C` such that for `J_E_dκ(ξ)` being the derivative w.r.t. `κ`
of
```
J_E(ξ) = -(1 + im * δ) / (1 - im * ϵ) * E(ξ) * inv(W(ξ))
```
we have
```
abs(J_E_dκ(ξ)) <= C * log(ξ) * ξ^(1 / λ.σ - 1)
```
for `ξ >= ξ₁`.
"""
function C_J_E_dκ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    f = ξ -> log(ξ) * ξ^(1 / λ.σ - 1)

    # FIXME: This is only an approximation. We multiply with 2 since
    # numerically it seems that the quotient is increasing in ξ.
    return 2 * abs(J_E(ξ₁, (λ, ArbSeries((κ, 1))))[1]) / f(ξ₁)
end
