_a(κ, ϵ, λ::CGLParams) = (1 / λ.σ + im * λ.ω / κ) / 2
function _a(κ::Arb, ϵ::Arb, λ::CGLParams{Arb})
    a = Acb()
    Arblib.inv!(Arblib.realref(a), λ.σ)
    Arblib.div!(Arblib.imagref(a), λ.ω, κ)
    return Arblib.mul_2exp!(a, a, -1)
end

_b(κ, ϵ, λ::CGLParams{T}) where {T} = Complex{T}(λ.d) / 2
_b(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}) = Acb(λ.d // 2)

_c(κ, ϵ, λ::CGLParams) = κ / 2(ϵ + im)
function _c(κ::Arb, ϵ::Arb, λ::CGLParams{Arb})
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

_abc(κ, ϵ, λ::CGLParams) = _a(κ, ϵ, λ), _b(κ, ϵ, λ), _c(κ, ϵ, λ)
_abc_dκ(κ, ϵ, λ::CGLParams) =
    _a(κ, ϵ, λ), _a_dκ(κ, ϵ, λ), _b(κ, ϵ, λ), _c(κ, ϵ, λ), _c_dκ(κ, ϵ, λ)
_abc_dϵ(κ, ϵ, λ::CGLParams) = _a(κ, ϵ, λ), _b(κ, ϵ, λ), _c(κ, ϵ, λ), _c_dϵ(κ, ϵ, λ)

function B_W(κ, ϵ, λ::CGLParams)
    (; δ) = λ
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(Arblib.imagref(Arblib.ref(c, 0)))
    elseif c isa Acb
        sign(Arblib.imagref(c))
    else
        sign(imag(c))
    end

    return _complex(-δ, 1) / κ * exp(-sgn * im * (b - a) * π) * c^b
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

function B_W_dϵ(κ, ϵ, λ::CGLParams)
    (; δ) = λ
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(Arblib.imagref(Arblib.ref(c, 0)))
    elseif c isa Acb
        sign(Arblib.imagref(c))
    else
        sign(imag(c))
    end

    return _complex(-δ, 1) / κ * exp(-sgn * im * (b - a) * π) * b * c_dϵ * c^(b - 1)
end

B_W_dϵ(κ, ϵ, λ) = ForwardDiff.derivative(ϵ -> B_W(κ, ϵ, λ), ϵ)
