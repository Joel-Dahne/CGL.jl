_a(κ, ϵ, λ::CGLParams) = (1 / λ.σ + im * λ.ω / κ) / 2
function _a(κ::Arb, ϵ::Arb, λ::CGLParams{Arb})
    # Acb(1 / λ.σ, λ.ω / κ) / 2
    a = Acb()
    Arblib.inv!(Arblib.realref(a), λ.σ)
    Arblib.div!(Arblib.imagref(a), λ.ω, κ)
    return Arblib.mul_2exp!(a, a, -1)
end

_b(κ, ϵ, λ::CGLParams{T}) where {T} = convert(Complex{T}, λ.d) / 2
_b(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}) = Acb(λ.d // 2)

_c(κ, ϵ, λ::CGLParams) = κ / 2(ϵ + im)
function _c(κ::Arb, ϵ::Arb, λ::CGLParams{Arb})
    # κ / 2Acb(ϵ, 1)
    c = Acb(κ)
    Arblib.mul_2exp!(c, c, -1)
    return Arblib.div!(c, c, Acb(ϵ, 1))
end

_a_dκ(κ, ϵ, λ::CGLParams) = -im * (λ.ω / κ^2) / 2
_a_dκ(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}) = Acb(0, -1) * (λ.ω / κ^2) / 2

_c_dκ(κ, ϵ, λ::CGLParams) = 1 / 2(ϵ + im)
_c_dκ(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}) = 1 / 2Acb(ϵ, 1)

_c_dϵ(κ, ϵ, λ::CGLParams) = -κ / 2(ϵ + im)^2
_c_dϵ(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}) = -κ / 2Acb(ϵ, 1)^2

function _abc(κ, ϵ, λ::CGLParams)
    a = _a(κ, ϵ, λ)
    b = _b(κ, ϵ, λ)
    c = _c(κ, ϵ, λ)

    return a, b, c
end

function _abc_dκ(κ, ϵ, λ::CGLParams{T}) where {T}
    a = _a(κ, ϵ, λ)
    a_dκ = _a_dκ(κ, ϵ, λ)
    b = _b(κ, ϵ, λ)
    c = _c(κ, ϵ, λ)
    c_dκ = _c_dκ(κ, ϵ, λ)

    return a, a_dκ, b, c, c_dκ
end

function _abc_dϵ(κ, ϵ, λ::CGLParams{T}) where {T}
    a = _a(κ, ϵ, λ)
    b = _b(κ, ϵ, λ)
    c = _c(κ, ϵ, λ)
    c_dϵ = _c_dϵ(κ, ϵ, λ)

    return a, b, c, c_dϵ
end

function B_W(κ, ϵ, λ::CGLParams{T}) where {T}
    (; δ) = λ

    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(Arblib.imagref(Arblib.ref(c, 0)))
    elseif c isa Acb
        sign(Arblib.imagref(c))
    else
        sign(imag(c))
    end

    if T == Arb
        return Acb(-δ, 1) / κ * exp(-sgn * im * (b - a) * π) * c^b
    else
        return (-δ + im) / κ * exp(-sgn * im * (b - a) * π) * c^b
    end
end

function B_W_dκ(κ::Arb, ϵ::Arb, λ::CGLParams)
    (; δ) = λ

    κ_series = ArbSeries((κ, 1))

    a, b, c = _abc(κ_series, ϵ, λ)

    sgn = sign(Arblib.imagref(Arblib.ref(c, 0)))

    res = Acb(-δ, 1) / κ_series * exp(-sgn * im * (b - a) * π) * c^b

    return res[1]
end

B_W_dκ(κ, ϵ, λ) = ForwardDiff.derivative(κ -> B_W(κ, ϵ, λ), κ)

function B_W_dϵ(κ, ϵ, λ::CGLParams{T}) where {T}
    (; δ) = λ

    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(Arblib.imagref(Arblib.ref(c, 0)))
    elseif c isa Acb
        sign(Arblib.imagref(c))
    else
        sign(imag(c))
    end

    if T == Arb
        return Acb(-δ, 1) / κ * exp(-sgn * im * (b - a) * π) * b * c_dϵ * c^(b - 1)
    else
        return (-δ + im) / κ * exp(-sgn * im * (b - a) * π) * b * c_dϵ * c^(b - 1)
    end
end
