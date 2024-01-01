# Coefficients in asymptotic expansion

function p_U!(res::Acb, k::Integer, a::Acb, b::Acb, z::Acb)
    tmp = a - b
    Arblib.add!(tmp, tmp, 1)
    Arblib.rising!(res, tmp, convert(UInt, k))
    Arblib.rising!(tmp, a, convert(UInt, k))
    Arblib.mul!(res, res, tmp)
    Arblib.neg!(tmp, z)
    Arblib.pow!(tmp, tmp, k)
    Arblib.mul!(tmp, tmp, factorial(k))
    return Arblib.div!(res, res, tmp)
end

p_U(k::Integer, a, b, z) = rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-z)^k)
p_U(k::Integer, a::Acb, b::Acb, z::Acb) = p_U!(zero(z), k, a, b, z)

function p_U_da!(res::Acb, k::Integer, a::Acb, b::Acb, z::Acb)
    tmp1 = zero(res)
    tmp2 = a - b
    Arblib.add!(tmp2, tmp2, 1)
    rising2!(res, tmp1, tmp2, convert(UInt, k), precision(res))

    tmp3 = zero(res)
    rising2!(tmp2, tmp3, a, convert(UInt, k), precision(res))

    Arblib.mul!(res, res, tmp3)
    Arblib.addmul!(res, tmp1, tmp2)

    # Reuse rising_1
    Arblib.neg!(tmp1, z)
    Arblib.pow!(tmp1, tmp1, k)
    Arblib.mul!(tmp1, tmp1, factorial(k))
    return Arblib.div!(res, res, tmp1)
end

p_U_da(k::Integer, a::Acb, b::Acb, z::Acb) = p_U_da!(zero(z), k, a, b, z)

# Bound for remainder terms in asymptotic expansion

function C_R_U(n::Integer, a, b, z)
    abs(imag(z)) > abs(imag(b - 2a)) ||
        throw(ArgumentError("assuming abs(imag(z)) > abs(imag(b - 2a))"))

    s = abs(b - 2a) / abs(z)
    ρ = abs(a^2 - a * b + b / 2) + s * (1 + s / 4) / (1 - s)^2

    # Bound for remainder for U
    # See https://fungrim.org/entry/461a54/
    return abs(rising(a, n) * rising(a - b + 1, n) / factorial(n)) *
           2sqrt(1 + Arb(π) * n / 2) / (1 - s) * exp(π * ρ / ((1 - s) * abs(z)))
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
