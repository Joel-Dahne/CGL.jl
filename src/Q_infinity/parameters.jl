_a(κ, ϵ, Λ::CGLParams) = (1 / Λ.σ + im * Λ.ω / κ) / 2
function _a(κ::Arb, ϵ::Arb, Λ::CGLParams{Arb})
    a = Acb()
    Arblib.inv!(Arblib.realref(a), Λ.σ)
    Arblib.div!(Arblib.imagref(a), Λ.ω, κ)
    return Arblib.mul_2exp!(a, a, -1)
end

_b(κ, ϵ, Λ::CGLParams{T}) where {T} = Complex{T}(Λ.d) / 2
_b(κ::Arb, ϵ::Arb, Λ::CGLParams{Arb}) = Acb(Λ.d // 2)

_c(κ, ϵ, Λ::CGLParams) = κ / 2(ϵ + im)
function _c(κ::Arb, ϵ::Arb, Λ::CGLParams{Arb})
    c = Acb(κ)
    Arblib.mul_2exp!(c, c, -1)
    return Arblib.div!(c, c, Acb(ϵ, 1))
end

_a_dκ(κ, ϵ, Λ::CGLParams) = -im * (Λ.ω / κ^2) / 2
_a_dκ(κ::Arb, ϵ::Arb, Λ::CGLParams{Arb}) = Acb(0, -1) * (Λ.ω / κ^2) / 2

_c_dκ(κ, ϵ, Λ::CGLParams) = 1 / 2(ϵ + im)
_c_dκ(κ::Arb, ϵ::Arb, Λ::CGLParams{Arb}) = 1 / 2Acb(ϵ, 1)

_c_dϵ(κ, ϵ, Λ::CGLParams) = -κ / 2(ϵ + im)^2
_c_dϵ(κ::Arb, ϵ::Arb, Λ::CGLParams{Arb}) = -κ / 2Acb(ϵ, 1)^2

_abc(κ, ϵ, Λ::CGLParams) = _a(κ, ϵ, Λ), _b(κ, ϵ, Λ), _c(κ, ϵ, Λ)
_abc_dκ(κ, ϵ, Λ::CGLParams) =
    _a(κ, ϵ, Λ), _a_dκ(κ, ϵ, Λ), _b(κ, ϵ, Λ), _c(κ, ϵ, Λ), _c_dκ(κ, ϵ, Λ)
_abc_dϵ(κ, ϵ, Λ::CGLParams) = _a(κ, ϵ, Λ), _b(κ, ϵ, Λ), _c(κ, ϵ, Λ), _c_dϵ(κ, ϵ, Λ)

function B_W(κ, ϵ, Λ::CGLParams)
    (; δ) = Λ
    a, b, c = _abc(κ, ϵ, Λ)

    sgn = if c isa AcbSeries
        sign(Arblib.imagref(Arblib.ref(c, 0)))
    elseif c isa Acb
        sign(Arblib.imagref(c))
    else
        sign(imag(c))
    end

    return _complex(-δ, 1) / κ * exp(-sgn * im * (b - a) * π) * c^b
end

function B_W_dκ(κ::Arb, ϵ::Arb, Λ::CGLParams)
    (; δ) = Λ
    κ_series = ArbSeries((κ, 1))
    a, b, c = _abc(κ_series, ϵ, Λ)

    sgn = sign(Arblib.imagref(Arblib.ref(c, 0)))

    res = Acb(-δ, 1) / κ_series * exp(-sgn * im * (b - a) * π) * c^b

    return res[1]
end

B_W_dκ(κ, ϵ, Λ) = ForwardDiff.derivative(κ -> B_W(κ, ϵ, Λ), κ)

function B_W_dϵ(κ::Arb, ϵ::Arb, Λ::CGLParams)
    (; δ) = Λ
    a, b, c, c_dϵ = _abc_dϵ(κ, ϵ, Λ)

    sgn = if c isa AcbSeries
        sign(Arblib.imagref(Arblib.ref(c, 0)))
    elseif c isa Acb
        sign(Arblib.imagref(c))
    else
        sign(imag(c))
    end

    return _complex(-δ, 1) / κ * exp(-sgn * im * (b - a) * π) * b * c_dϵ * c^(b - 1)
end

B_W_dϵ(κ, ϵ, Λ) = ForwardDiff.derivative(ϵ -> B_W(κ, ϵ, Λ), ϵ)
