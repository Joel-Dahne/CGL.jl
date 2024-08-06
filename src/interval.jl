# Conversions between Arb and BareInterval/Interval

Arblib.Arb(x::Union{Interval,BareInterval}) =
    isempty_interval(x) ? indeterminate(Arb) : Arb((inf(x), sup(x)))

IntervalArithmetic.bareinterval(::Type{T}, x::Arb) where {T} =
    isnan(x) ? nai(T).bareinterval : bareinterval(T, getinterval(BigFloat, x)...)

IntervalArithmetic.interval(::Type{T}, x::Arb) where {T} =
    isnan(x) ? nai(T) : interval(T, getinterval(BigFloat, x)...)

IntervalArithmetic.interval(::Type{T}, (x, y)::NTuple{2,Arb}) where {T} =
    (isnan(Arblib.midref(x)) || isnan(Arblib.midref(y))) ? nai(T) :
    interval(T, BigFloat(lbound(x)), BigFloat(ubound(y)))

IntervalArithmetic.bareinterval(x::Arb) = bareinterval(Float64, x)
IntervalArithmetic.interval(x::Arb) = interval(Float64, x)
IntervalArithmetic.interval((x, y)::NTuple{2,Arb}) = interval(Float64, (x, y))

Base.convert(::Type{Arb}, x::BareInterval) = Arb(x)
Base.convert(::Type{BareInterval{T}}, x::Arb) where {T} = bareinterval(T, x)
