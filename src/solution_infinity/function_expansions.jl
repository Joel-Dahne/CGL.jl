# Coefficients in asymptotic expansion

p_U(k::Integer, a, b) = rising(a, k) * rising(a - b + 1, k) / factorial(k)

p_U_da(k::Integer, a::Acb, b::Acb) =
    let a = AcbSeries((a, 1))
        res = rising(a, k) * rising(a - b + 1, k) / factorial(k)
        return res[1]
    end

# Bound for remainder terms in asymptotic expansion

function C_R_U(n::Integer, a, b, z)
    abs(imag(z)) > abs(imag(b - 2a)) ||
        throw(ArgumentError("assuming abs(imag(z)) > abs(imag(b - 2a))"))

    s = abs(b - 2a) / abs(z)
    ρ = abs(a^2 - a * b + b / 2) + s * (1 + s / 4) / (1 - s)^2

    # Bound for remainder for U
    # See https://fungrim.org/entry/461a54/
    return abs(p_U(n, a, b)) * 2sqrt(1 + Arb(π) * n / 2) / (1 - s) *
           exp(π * ρ / ((1 - s) * abs(z)))
end

function C_R_U_1(n::Integer, a, b, z)
    @assert 0 < real(a) < real(b)
    @assert 0 < real(a - b + n + 1)

    γ = Arblib.arg!(Arb(), z)

    ρ_γ = if real(z) >= 0 # Corresponds to abs(γ) <= π / 2
        Arb(1)
    else
        # This is always >= 1, so correct even when z overlaps the
        # imaginary axis.
        inv(sin(π - abs(γ)))
    end

    C_χ =
        ρ_γ * abs(sinpi(a - b + 1)) / π * gamma(-real(a - b)) * gamma(real(a - b + n + 1)) /
        factorial(n)

    C1 = gamma(real(a + n))

    C2 = inv(real(a + n)^2) + gamma(real(a + n + 1))

    return C_χ / abs(gamma(a)) * exp(π * imag(n + a)) / (1 + π / abs(log(abs(z)))) *
           (C1 + (abs(γ) * C1 + C2) / abs(log(abs(z))))
end

function C_R_U_2(n::Integer, a, b, z)
    @assert 0 < real(a) < real(b)
    @assert 0 < real(a - b + n + 1)

    γ = Arblib.arg!(Arb(), z)

    ρ_γ = if real(z) >= 0 # Corresponds to abs(γ) <= π / 2
        Arb(1)
    else
        # This is always >= 1, so correct even when z overlaps the
        # imaginary axis.
        inv(sin(π - abs(γ)))
    end

    C_χ_a =
        ρ_γ * abs(cospi(a - b + 1)) * gamma(-real(a - b)) * gamma(real(a - b + n + 1)) /
        factorial(n) +
        ρ_γ * abs(sinpi(a - b + 1)) / π * (
            inv(real(a - b)^2) +
            gamma(-real(a - b)) * gamma(real(a - b + n)) / factorial(n - 1)
        )

    C1 = gamma(real(a + n))

    return C_χ_a * C1 / abs(gamma(a)) * exp(π * imag(n + a)) / abs(log(abs(z)))
end
