# Coefficients in asymptotic expansion

p_U(k::Integer, a, b) = rising(a, k) * rising(a - b + 1, k) / factorial(k)

p_U_da(k::Integer, a::Acb, b::Acb) =
    let a = AcbSeries((a, 1))
        res = rising(a, k) * rising(a - b + 1, k) / factorial(k)
        return res[1]
    end

function p_P(k::Integer, κ, λ::CGLParams)
    a, b, c = _abc(κ, λ)

    return c^-a * rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k)
end

function p_E(k::Integer, κ, λ::CGLParams)
    a, b, c = _abc(κ, λ)

    return (-c)^(a - b) * rising(b - a, k) * rising(b - 2a + 1, k) / (factorial(k) * c^k)
end

# Bound for remainder terms in asymptotic expansion

function C_R_U(n::Integer, a, b, z)
    s = abs(b - 2a) / abs(z)
    ρ = abs(a^2 - a * b + b / 2) + s * (1 + s / 4) / (1 - s)^2

    # Bound for remainder for hypgeom_u
    # See https://fungrim.org/entry/461a54/
    return abs(p_U(n, a, b)) * 2sqrt(1 + Arb(π) * n / 2) / (1 - s) *
           exp(π * ρ / ((1 - s) * abs(z)))
end

# TODO: Implement proper bound. Now we take 10 times the magnitude of
# the corresponding term
function C_R_U_da(n::Integer, a, b, z)
    return 10abs(p_U_da(n, a, b))
end
