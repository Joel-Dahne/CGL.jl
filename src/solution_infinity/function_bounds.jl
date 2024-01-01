"""
    C_U(a::Acb, b::Acb, z₁::Acb)

Return `C` such that
```
abs(U(a, b, z)) <= C * abs(z^-a)
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
function C_U(a::Acb, b::Acb, z₁::Acb, n::Integer = 5)
    S = sum(0:n-1) do k
        abs(p_U(k, a, b, z₁))
    end

    return S + C_R_U(n, a, b, z₁) * abs(z₁)^-n
end

C_U_dz(a::Acb, b::Acb, z₁::Acb, n::Integer = 1) =
    if n < 0
        throw(ArgumentError("n must be non-negative"))
    elseif n == 0
        return C_U(a, b, z₁)
    else
        return C_U(a + n, b + n, z₁) * abs(rising(a, n))
    end

function C_U_da(a::Acb, b::Acb, z₁::Acb, n::Integer = 5)
    S1 = sum(0:n-1) do k
        abs(p_U(k, a, b, z₁))
    end

    S2 = sum(0:n-1) do k
        abs(p_U_da(k, a, b, z₁))
    end

    # Note that Γ'(a) / Γ(a) is exactly the digamma function
    R =
        (1 + abs(digamma(a) / log(z₁))) * C_R_U(n, a, b, z₁) +
        C_R_U_1(n, a, b, z₁) +
        C_R_U_2(n, a, b, z₁)

    return S1 + S2 / abs(log(z₁)) + R * abs(z₁)^-n
end

function C_P(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    return C_U(a, b, c * ξ₁^2) * abs(c^-a)
end

function C_E(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    return C_U(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b))
end

function C_P_dξ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    return abs(2c^-a) * C_U_dz(a, b, c * ξ₁^2)
end

function C_P_dξ_dξ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)
    z₁ = c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs(2(a + 1) * p_U(k, a + 2, b + 2, z₁) - p_U(k, a + 1, b + 1, z₁))
    end

    R = 2abs(a + 1) * C_R_U(n, a + 2, b + 2, z₁) + C_R_U(n, a + 1, b + 1, z₁)

    return abs(2a * c^-a) * (S + R * abs(z₁)^-n)
end

function C_P_dξ_dξ_dξ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)
    z₁ = c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs(-2(a + 2) * p_U(k, a + 3, b + 3, z₁) + 3p_U(k, a + 2, b + 2, z₁))
    end

    R = 2abs(a + 2) * C_R_U(n, a + 3, b + 3, z₁) + 3C_R_U(n, a + 2, b + 2, z₁)

    return abs(4a * (a + 1) * c^-a) * (S + R * abs(z₁)^-n)
end

function C_E_dξ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    C1 = abs((-c)^(a - b)) * C_U(b - a, b, -c * ξ₁^2)
    C2 = abs((-c)^(a - b - 1)) * C_U_dz(b - a, b, -c * ξ₁^2)

    return abs(2c) * C1 + abs(2c) * C2 * ξ₁^-2
end

function C_E_dξ_dξ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    C1 = abs((-c)^(a - b)) * C_U(b - a, b, -c * ξ₁^2)
    C2 = abs((-c)^(a - b - 1)) * C_U_dz(b - a, b, -c * ξ₁^2)
    C3 = abs((-c)^(a - b - 2)) * C_U_dz(b - a, b, -c * ξ₁^2, 2)

    return abs(4c^2 + 2c * ξ₁^-2) * C1 +
           abs(8c^2 + 2c * ξ₁^-2) * C2 * ξ₁^-2 +
           abs(4c^2) * C3 * ξ₁^-4
end

function C_E_dξ_dξ_dξ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    C1 = abs((-c)^(a - b)) * C_U(b - a, b, -c * ξ₁^2)
    C2 = abs((-c)^(a - b - 1)) * C_U_dz(b - a, b, -c * ξ₁^2)
    C3 = abs((-c)^(a - b - 2)) * C_U_dz(b - a, b, -c * ξ₁^2, 2)
    C4 = abs((-c)^(a - b - 3)) * C_U_dz(b - a, b, -c * ξ₁^2, 3)

    return abs(8c^3 + 12c^2 * ξ₁^-2) * C1 +
           abs(24c^3 + 24c^2 * ξ₁^-2) * C2 * ξ₁^-2 +
           abs(24c^3 + 12c^2 * ξ₁^-2) * C3 * ξ₁^-4 +
           abs(8c^3) * C4 * ξ₁^-6
end

# IMPROVE: This upper bound can be improved
function C_P_dκ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    C1 = C_U_dz(a, b, c * ξ₁^2) * abs(c^(-a - 1) * c_dκ)

    C2 = C_U_da(a, b, c * ξ₁^2) * abs(c^-a * a_dκ) * (2 + abs(log(c)) / log(ξ₁))

    return C1 / log(ξ₁) + C2
end

function C_E_dκ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    C1 = C_U(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b) * c_dκ)

    C2 =
        C_U_da(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b) * a_dκ) * (2 + abs(log(c)) / log(ξ₁))

    C3 = C_U_dz(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b - 1) * c_dκ)

    return C1 + C2 * log(ξ₁) * ξ₁^-2 + C3 * ξ₁^-2
end

# IMPROVE: This upper bound can be improved
function C_P_dξ_dκ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    C1 = C_U_dz(a, b, c * ξ₁^2) * abs(c^(-a - 1) * 2c_dκ)

    C2 = C_U_dz(a, b, c * ξ₁^2, 2) * abs(c^(-a - 1) * 2c_dκ)

    C3 =
        C_U_da(a + 1, b + 1, c * ξ₁^2) * abs(c^-a * 2a_dκ * a) * (2 + abs(log(c)) / log(ξ₁))

    C4 = C_U(a + 1, b + 1, c * ξ₁^2) * abs(c^-a * 2a_dκ)

    return C1 / log(ξ₁) + C2 / log(ξ₁) + C3 + C4 / log(ξ₁)
end

function C_E_dξ_dκ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)
    (; d) = λ

    C1 = 2C_U(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b) * c_dκ) * (abs(c) + ξ₁^-2)

    C2 = 2C_U_dz(b - a, b, -c * ξ₁^2) * abs((-c)^(a - b - 1) * c_dκ) * (abs(2c) + ξ₁^-2)

    C3 =
        2C_U_da(b - a, b, -c * ξ₁^2) *
        abs((-c)^(a - b) * a_dκ * c) *
        (2 + abs(log(-c)) / log(ξ₁))

    C4 = 2C_U_dz(b - a, b, -c * ξ₁^2, 2) * abs((-c)^(a - b - 2) * c * c_dκ)

    C5 = 2C_U(b - a + 1, b + 1, -c * ξ₁^2) * abs((-c)^(a - b - 1) * a_dκ * c)

    C6 =
        2C_U_da(b - a + 1, b + 1, -c * ξ₁^2) *
        abs((b - a) * (-c)^(a - b - 1) * a_dκ * c) *
        (2 + abs(log(-c)) / log(ξ₁))

    return C1 +
           C2 * ξ₁^-2 +
           C3 * log(ξ₁) * ξ₁^-2 +
           C4 * ξ₁^-2 +
           C5 * ξ₁^-4 +
           C6 * log(ξ₁) * ξ₁^-4
end

# IMPROVE: This upper bound can be improved
function C_P_dξ_dξ_dκ(κ::Arb, λ::CGLParams{Arb}, ξ₁::Arb)
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    C1 = C_U_dz(a, b, c * ξ₁^2) * abs(c^(-a - 1) * 2c_dκ)

    C2 = C_U_dz(a, b, c * ξ₁^2, 2) * abs(c^(-a - 1) * 10c_dκ)

    C3 = C_U_dz(a, b, c * ξ₁^2, 3) * abs(c^(-a - 1) * 4c_dκ)

    C4 = C_U(a + 1, b + 1, c * ξ₁^2) * abs(c^-a * 2a_dκ)

    C5 =
        C_U_da(a + 1, b + 1, c * ξ₁^2) * abs(c^-a * 2a_dκ * a) * (2 + abs(log(c)) / log(ξ₁))

    C6 = C_U(a + 2, b + 2, c * ξ₁^2) * abs(c^-a * 4a_dκ * (2a + 1))

    C7 =
        C_U_da(a + 2, b + 2, c * ξ₁^2) *
        abs(c^-a * 4a_dκ * a * (a + 1)) *
        (2 + abs(log(c)) / log(ξ₁))

    return C1 / log(ξ₁) +
           C2 / log(ξ₁) +
           C3 / log(ξ₁) +
           C4 / log(ξ₁) +
           C5 +
           C6 / log(ξ₁) +
           C7
end

C_J_P(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb}) = abs(B_W(κ, λ)) * C_P(κ, λ, ξ₁)

C_J_E(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb}) = abs(B_W(κ, λ)) * C_E(κ, λ, ξ₁)

function C_J_P_dξ(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    _, _, c = _abc(κ, λ)
    (; d) = λ

    return abs(B_W(κ, λ)) *
           (C_P(κ, λ, ξ₁) * (abs(2c) + (d - 1) * ξ₁^-2) + C_P_dξ(κ, λ, ξ₁) * ξ₁^-2)
end

function C_J_E_dξ(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, λ)
    (; d) = λ
    z₁ = -c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs((d - 1) * p_U(k, b - a, b, z₁) - 2(b - a) * p_U(k, b - a + 1, b + 1, z₁))
    end

    R = (d - 1) * C_R_U(n, b - a, b, z₁) + abs(2(b - a)) * C_R_U(n, b - a + 1, b + 1, z₁)

    return abs(B_W(κ, λ) * (-c)^(a - b)) * (S + R * abs(z₁)^-n)
end

function C_J_P_dξ_dξ(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    _, _, c = _abc(κ, λ)
    (; d) = λ

    return abs(B_W(κ, λ)) * (
        C_P(κ, λ, ξ₁) *
        (abs(4c^2) + abs(2c) * (2d - 1) * ξ₁^-2 + (d - 1) * (d - 2) * ξ₁^-4) +
        C_P_dξ(κ, λ, ξ₁) * (abs(4c) + 2(d - 1) * ξ₁^-2) * ξ₁^-2 +
        C_P_dξ_dξ(κ, λ, ξ₁) * ξ₁^-4
    )
end

function C_J_E_dξ_dξ(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, λ)
    (; d) = λ
    z₁ = -c * ξ₁^2
    n = 5

    S = sum(0:n-1) do k
        abs(
            (d - 1) * (d - 2) * p_U(k, b - a, b, z₁) -
            2(2d - 1) * (b - a) * p_U(k, b - a + 1, b + 1, z₁) +
            4(b - a) * (b - a + 1) * p_U(k, b - a + 2, b + 2, z₁),
        )
    end

    R =
        (d - 1) * (d - 2) * C_R_U(n, b - a, b, z₁) +
        abs(2(2d - 1) * (b - a)) * C_R_U(n, b - a + 1, b + 1, z₁) +
        abs(4(b - a) * (b - a + 1)) * C_R_U(n, b - a + 2, b + 2, z₁)

    return abs(B_W(κ, λ) * (-c)^(a - b)) * (S + R * abs(z₁)^-n)
end

function C_J_P_dκ(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)
    (; d) = λ

    C1 = C_P(κ, λ, ξ₁) * (abs(c_dκ * B_W(κ, λ)) + abs(B_W_dκ(κ, λ) * ξ₁^-2))

    C2 = C_P_dκ(κ, λ, ξ₁) * abs(B_W(κ, λ))

    return C1 + C2 * log(ξ₁) * ξ₁^-2
end

# IMPROVE: This upper bound can be improved
function C_J_E_dκ(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)
    (; d) = λ

    C1 = abs(B_W_dκ(κ, λ) * (-c)^(a - b)) * C_U(b - a, b, -c * ξ₁^2)

    C2 =
        abs(B_W(κ, λ) * (-c)^(a - b) * a_dκ) *
        (2 + abs(log(-c)) / log(ξ₁)) *
        C_U_da(b - a, b, -c * ξ₁^2)

    C3 = abs(B_W(κ, λ) * (-c)^(a - b - 1) * c_dκ) * C_U_dz(b - a, b, -c * ξ₁^2)

    return C1 / log(ξ₁) + C2 + C1 / log(ξ₁)
end

function C_D(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    C1 = abs(c_dκ * B_W(κ, λ)) * C_P(κ, λ, ξ₁)

    C2 = abs(B_W_dκ(κ, λ)) * C_P(κ, λ, ξ₁)

    C3 = abs(B_W(κ, λ)) * C_P_dκ(κ, λ, ξ₁)

    return C1 + (C2 + C3 * log(ξ₁)) * ξ₁^-2
end

function C_D_dξ(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    C1 = abs(c_dκ * B_W(κ, λ)) * C_P_dξ(κ, λ, ξ₁)

    C2 = abs(B_W_dκ(κ, λ)) * C_P_dξ(κ, λ, ξ₁)

    C3 = abs(2B_W_dκ(κ, λ)) * C_P(κ, λ, ξ₁)

    C4 = abs(B_W(κ, λ)) * C_P_dξ_dκ(κ, λ, ξ₁)

    C5 = abs(2B_W(κ, λ)) * C_P_dκ(κ, λ, ξ₁)

    return C1 + (C2 + C3 + (C4 + C5) * log(ξ₁)) * ξ₁^-2
end

function C_D_dξ_dξ(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, a_dκ, b, c, c_dκ = _abc_dκ(κ, λ)

    C1 = abs(c_dκ * B_W(κ, λ)) * C_P_dξ_dξ(κ, λ, ξ₁)

    C2 = abs(B_W_dκ(κ, λ)) * C_P_dξ_dξ(κ, λ, ξ₁)

    C3 = abs(4B_W_dκ(κ, λ)) * C_P_dξ(κ, λ, ξ₁)

    C4 = abs(6B_W_dκ(κ, λ)) * C_P(κ, λ, ξ₁)

    C5 = abs(B_W(κ, λ)) * C_P_dξ_dξ_dκ(κ, λ, ξ₁)

    C6 = abs(2B_W(κ, λ)) * C_P_dξ_dκ(κ, λ, ξ₁)

    C7 = abs(6B_W(κ, λ)) * C_P_dκ(κ, λ, ξ₁)

    return C1 + (C2 + C3 + C4 + (C5 + C6 + C7) * log(ξ₁)) * ξ₁^-2
end

struct FunctionBounds
    P::Arb
    P_dξ::Arb
    P_dξ_dξ::Arb
    P_dξ_dξ_dξ::Arb
    P_dκ::Arb
    P_dξ_dκ::Arb
    P_dξ_dξ_dκ::Arb
    E::Arb
    E_dξ::Arb
    E_dξ_dξ::Arb
    E_dξ_dξ_dξ::Arb
    E_dκ::Arb
    E_dξ_dκ::Arb
    J_P::Arb
    J_P_dξ::Arb
    J_P_dξ_dξ::Arb
    J_P_dκ::Arb
    J_E::Arb
    J_E_dξ::Arb
    J_E_dξ_dξ::Arb
    J_E_dκ::Arb
    D::Arb
    D_dξ::Arb
    D_dξ_dξ::Arb

    function FunctionBounds(κ::Arb, ξ₁::Arb, λ::CGLParams{Arb}; skip_dκ::Bool = false)
        return new(
            C_P(κ, λ, ξ₁),
            C_P_dξ(κ, λ, ξ₁),
            C_P_dξ_dξ(κ, λ, ξ₁),
            C_P_dξ_dξ_dξ(κ, λ, ξ₁),
            skip_dκ ? indeterminate(κ) : C_P_dκ(κ, λ, ξ₁),
            skip_dκ ? indeterminate(κ) : C_P_dξ_dκ(κ, λ, ξ₁),
            skip_dκ ? indeterminate(κ) : C_P_dξ_dξ_dκ(κ, λ, ξ₁),
            C_E(κ, λ, ξ₁),
            C_E_dξ(κ, λ, ξ₁),
            C_E_dξ_dξ(κ, λ, ξ₁),
            C_E_dξ_dξ_dξ(κ, λ, ξ₁),
            skip_dκ ? indeterminate(κ) : C_E_dκ(κ, λ, ξ₁),
            skip_dκ ? indeterminate(κ) : C_E_dξ_dκ(κ, λ, ξ₁),
            C_J_P(κ, ξ₁, λ),
            C_J_P_dξ(κ, ξ₁, λ),
            C_J_P_dξ_dξ(κ, ξ₁, λ),
            skip_dκ ? indeterminate(κ) : C_J_P_dκ(κ, ξ₁, λ),
            C_J_E(κ, ξ₁, λ),
            C_J_E_dξ(κ, ξ₁, λ),
            C_J_E_dξ_dξ(κ, ξ₁, λ),
            skip_dκ ? indeterminate(κ) : C_J_E_dκ(κ, ξ₁, λ),
            skip_dκ ? indeterminate(κ) : C_D(κ, ξ₁, λ),
            skip_dκ ? indeterminate(κ) : C_D_dξ(κ, ξ₁, λ),
            skip_dκ ? indeterminate(κ) : C_D_dξ_dξ(κ, ξ₁, λ),
        )
    end
end
