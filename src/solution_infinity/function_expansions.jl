# Coefficient in asymptotic expansion of P
function p_P(k::Integer, κ, λ::AbstractGLParams)
    a, b, c = _abc(κ, λ)

    return c^-a * rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k)
end

# Bound for remainder in asymptotic expansion with k terms
# TODO: Better name for thus function
function C_P(k::Integer, κ::Arb, λ::AbstractGLParams{Arb}, ξ₁::Arb)
    a, b, c = _abc(κ, λ)

    z₁ = c * ξ₁^2
    s = abs(b - 2a) / abs(z₁)
    ρ = abs(a^2 - a * b + b / 2) + s * (1 + s / 4) / (1 - s)^2

    # Bound for remainder for hypgeom_u
    C_hypgeom_u =
        abs(rising(a, k) * rising(a - b + 1, k) / (factorial(k) * abs(z₁)^k)) *
        2sqrt(1 + Arb(π) * k / 2) / (1 - s) * exp(π * ρ / ((1 - s) * abs(z₁)))

    return C_hypgeom_u * abs(c^-a)
end
