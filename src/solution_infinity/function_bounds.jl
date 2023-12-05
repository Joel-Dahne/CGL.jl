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

# Extended help
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

    S = sum(0:n-1) do k
        abs(p_U(k, a, b) * z₁^-k)
    end

    return S + C_R_U(n, a, b, z₁) * abs(z₁)^-n
end

C_hypgeom_u_dz(a::Acb, b::Acb, z₁::Acb, n::Integer = 1) =
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return C_hypgeom_u(a, b, z₁)
    else
        return C_hypgeom_u(a + n, b + n, z₁) * abs(rising(a, n))
    end

function C_hypgeom_u_da(a::Acb, b::Acb, z₁::Acb, n::Integer = 5)
    abs(imag(z₁)) > abs(imag(b - 2a)) ||
        throw(ArgumentError("assuming abs(imag(z₁)) > abs(imag(b - 2a))"))

    C1 = C_hypgeom_u(a, b, z₁, n)

    C2 = let
        S = sum(0:n-1) do k
            abs(p_U_da(k, a, b) * (-z₁)^-k)
        end

        R = C_R_U_da(n, a, b, z₁)

        S + R * abs(z₁)^-n
    end

    return C1 + C2 / abs(log(z₁))
end

function C_P(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    return C_hypgeom_u(a, b, c * ξ₁^2) * abs(c^-a)
end

function C_E(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    return C_hypgeom_u(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b))
end

function C_P_dξ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    return abs(2c^-a) * C_hypgeom_u_dz(a, b, c * ξ₁^2)
end

function C_P_dξ_dξ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)
    z₁ = c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs((2(a + 1) * p_U(k, a + 2, b + 2) - p_U(k, a + 1, b + 1)) * z₁^-k)
    end

    R = 2abs(a + 1) * C_R_U(n, a + 2, b + 2, z₁) + C_R_U(n, a + 1, b + 1, z₁)

    return abs(2a * c^-a) * (S + R * abs(z₁)^-n)
end

function C_P_dξ_dξ_dξ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)
    z₁ = c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs((-2(a + 2) * p_U(k, a + 3, b + 3) + 3p_U(k, a + 2, b + 2)) * z₁^-k)
    end

    R = 2abs(a + 2) * C_R_U(n, a + 3, b + 3, z₁) + 3C_R_U(n, a + 2, b + 2, z₁)

    return abs(4a * (a + 1) * c^-a) * (S + R * abs(z₁)^-n)
end

function C_E_dξ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    C1 = abs((-c)^(a - b)) * C_hypgeom_u(b - a, b, -c * ξ₁^2)
    C2 = abs((-c)^(a - b - 1)) * C_hypgeom_u_dz(b - a, b, -c * ξ₁^2)

    return abs(2c) * C1 + abs(2c) * C2 * ξ₁^-2
end

function C_E_dξ_dξ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    C1 = abs((-c)^(a - b)) * C_hypgeom_u(b - a, b, -c * ξ₁^2)
    C2 = abs((-c)^(a - b - 1)) * C_hypgeom_u_dz(b - a, b, -c * ξ₁^2)
    C3 = abs((-c)^(a - b - 2)) * C_hypgeom_u_dz(b - a, b, -c * ξ₁^2, 2)

    return abs(4c^2 + 2c * ξ₁^-2) * C1 +
           abs(8c^2 + 2c * ξ₁^-2) * C2 * ξ₁^-2 +
           abs(4c^2) * C3 * ξ₁^-4
end

function C_E_dξ_dξ_dξ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    C1 = abs((-c)^(a - b)) * C_hypgeom_u(b - a, b, -c * ξ₁^2)
    C2 = abs((-c)^(a - b - 1)) * C_hypgeom_u_dz(b - a, b, -c * ξ₁^2)
    C3 = abs((-c)^(a - b - 2)) * C_hypgeom_u_dz(b - a, b, -c * ξ₁^2, 2)
    C4 = abs((-c)^(a - b - 3)) * C_hypgeom_u_dz(b - a, b, -c * ξ₁^2, 3)

    return abs(8c^3 + 12c^2 * ξ₁^-2) * C1 +
           abs(24c^3 + 24c^2 * ξ₁^-2) * C2 * ξ₁^-2 +
           abs(24c^3 + 12c^2 * ξ₁^-2) * C3 * ξ₁^-4 +
           abs(8c^3) * C4 * ξ₁^-6
end

# IMPROVE: This upper bound can be improved
function C_P_dκ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    C1 = C_hypgeom_u_da(a, b, c * ξ₁^2) * abs(2 + log(c) / log(ξ₁)) * abs(c^-a * a_dκ)

    C2 = C_hypgeom_u_dz(a, b, c * ξ₁^2) * abs(c^(-a - 1) * c_dκ)

    C = C1 + C2 / log(ξ₁)

    return C
end

function C_E_dκ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    C1 = C_hypgeom_u(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b) * c_dκ)

    C2 =
        C_hypgeom_u_da(b - a, b, -c * ξ₁^2) *
        abs(2 + log(c) / log(ξ₁)) *
        abs((-c)^(a - b) * a_dκ)

    C3 = C_hypgeom_u_dz(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b - 1) * c_dκ)

    C = C1 + C2 * log(ξ₁) * ξ₁^-2 + C3 * ξ₁^-2

    return C
end

# IMPROVE: This upper bound can be improved
function C_P_dξ_dκ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    C1 = C_hypgeom_u_dz(a, b, c * ξ₁^2) * abs(c^(-a - 1) * 2c_dκ)

    C2 = C_hypgeom_u(a + 1, b + 1, c * ξ₁^2) * abs(2c * c^(-a - 1) * a_dκ)

    C3 =
        abs(a) *
        C_hypgeom_u_da(a + 1, b + 1, c * ξ₁^2) *
        abs(2 + log(-c) / log(ξ₁)) *
        abs(2c * c^(-a - 1) * a_dκ)

    C4 = C_hypgeom_u_dz(a, b, c * ξ₁^2, 2) * abs(c^(-a - 2) * 2c * c_dκ)

    return C1 / log(ξ₁) + C2 / log(ξ₁) + C3 + C4 / log(ξ₁)
end

function C_E_dξ_dκ(κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)
    (; d) = λ

    C1 = 2C_hypgeom_u(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b) * c_dκ) * abs(c + ξ₁^-2)

    C2 =
        2C_hypgeom_u_dz(b - a, b, -c * ξ₁^2) *
        abs((-c)^(a - b - 1) * c_dκ) *
        abs(2c + ξ₁^-2)

    C3 =
        2C_hypgeom_u_da(b - a, b, -c * ξ₁^2) *
        abs(2 + log(-c) / log(ξ₁)) *
        abs((-c)^(a - b) * a_dκ * c)

    C4 = 2C_hypgeom_u_dz(b - a, b, -c * ξ₁^2, 2) * abs((-c)^(a - b - 2) * c * c_dκ)

    C5 = 2C_hypgeom_u(b - a + 1, b + 1, -c * ξ₁^2) * abs((-c)^(a - b - 1) * a_dκ * c)

    C6 =
        2C_hypgeom_u_da(b - a + 1, b + 1, -c * ξ₁^2) *
        abs(2 + log(-c) / log(ξ₁)) *
        abs((b - a) * (-c)^(a - b - 1) * a_dκ * c)

    return C1 +
           C2 * ξ₁^-2 +
           C3 * log(ξ₁) * ξ₁^-2 +
           C4 * ξ₁^-2 +
           C5 * ξ₁^-4 +
           C6 * log(ξ₁) * ξ₁^-4
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

C_J_P(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb}) = abs(B_W(κ, λ)) * C_P(κ, λ, ξ₁)

C_J_E(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb}) = abs(B_W(κ, λ)) * C_E(κ, λ, ξ₁)

function C_J_P_dξ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    _, _, c = _abc(κ, λ)
    (; d) = λ

    return abs(B_W(κ, λ)) *
           (C_P(κ, λ, ξ₁) * (abs(2c) + (d - 1) * ξ₁^-2) + C_P_dξ(κ, λ, ξ₁) * ξ₁^-2)
end

function C_J_E_dξ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    a, b, c = _abc(κ, λ)
    (; d) = λ
    z₁ = -c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs(((d - 1) * p_U(k, b - a, b) - 2(b - a) * p_U(k, b - a + 1, b + 1)) * (-z₁)^-k)
    end

    R = (d - 1) * C_R_U(n, b - a, b, z₁) + abs(2(b - a)) * C_R_U(n, b - a + 1, b + 1, z₁)

    return abs(B_W(κ, λ) * (-c)^(a - b)) * (S + R * abs(z₁)^-n)
end

function C_J_P_dξ_dξ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    _, _, c = _abc(κ, λ)
    (; d) = λ

    return abs(B_W(κ, λ)) * (
        C_P(κ, λ, ξ₁) *
        (abs(4c^2) + abs(2c) * (2d - 1) * ξ₁^-2 + (d - 1) * (d - 2) * ξ₁^-4) +
        C_P_dξ(κ, λ, ξ₁) * (abs(4c) + 2(d - 1) * ξ₁^-2) * ξ₁^-2 +
        C_P_dξ_dξ(κ, λ, ξ₁) * ξ₁^-4
    )
end

function C_J_E_dξ_dξ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    a, b, c = _abc(κ, λ)
    (; d) = λ
    z₁ = -c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs(
            (
                (d - 1) * (d - 2) * p_U(k, b - a, b) -
                2(2d - 1) * (b - a) * p_U(k, b - a + 1, b + 1) +
                4(b - a) * (b - a + 1) * p_U(k, b - a + 2, b + 2)
            ) * (-z₁)^-k,
        )
    end

    R =
        (d - 1) * (d - 2) * C_R_U(n, b - a, b, z₁) +
        abs(2(2d - 1) * (b - a)) * C_R_U(n, b - a + 1, b + 1, z₁) +
        abs(4(b - a) * (b - a + 1)) * C_R_U(n, b - a + 2, b + 2, z₁)

    return abs(B_W(κ, λ) * (-c)^(a - b)) * (S + R * abs(z₁)^-n)
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
    return 1.01 * abs(J_P_dκ(ξ₁, (λ, κ))) / f(ξ₁)
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
    return 2 * abs(J_E_dκ(ξ₁, (λ, κ))) / f(ξ₁)
end

function C_D(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    f = ξ -> ξ^(-1 / λ.σ)

    # FIXME: This is only an approximation.
    return 1.01 * abs(D(ξ₁, (λ, κ))) / f(ξ₁)
end

function C_D_dξ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    f = ξ -> ξ^(-1 / λ.σ - 1)

    # FIXME: This is only an approximation.
    return 1.01 * abs(D_dξ(ξ₁, (λ, κ))) / f(ξ₁)
end

function C_D_dξ_dξ(κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    f = ξ -> ξ^(-1 / λ.σ - 2)

    # FIXME: This is only an approximation.
    return 1.01 * abs(D_dξ_dξ(ξ₁, (λ, κ))) / f(ξ₁)
end
