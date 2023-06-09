# There is no version of this in Arblib
function Arblib.indeterminate!(x::Union{ArbSeries,AcbSeries})
    for i = 0:Arblib.degree(x)
        Arblib.indeterminate!(Arblib.ref(x, i))
    end
    # Since we manually set the coefficients of the polynomial we
    # need to also manually set the degree.
    Arblib.cstruct(x).length = Arblib.degree(x) + 1
    return x
end

Arblib.midpoint(::Type{Acb}, x::Arblib.AcbOrRef) =
    Acb(Arblib.midref(Arblib.realref(x)), Arblib.midref(Arblib.imagref(x)))

"""
    indeterminate(x)

Construct an indeterminate version of `x`.
"""
indeterminate(x::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}) = Arblib.indeterminate!(zero(x))
indeterminate(::Type{T}) where {T<:Union{Arb,Acb}} = Arblib.indeterminate!(zero(T))
indeterminate(x::Union{ArbSeries,AcbSeries}) = Arblib.indeterminate!(zero(x))

"""
    iswide(x; cutoff = 10)

Return true if `x` is wide in the meaning that the effective relative
accuracy of `x` measured in bits is more than `cutoff` lower than it's
precision. For `x` not of type `Arb` or `Acb` this always return
`false`. For `x` of type `ArbSeries` or `AcbSeries` it checks the
first coefficient.
"""
iswide(x::Union{Arblib.ArbOrRef,Arblib.AcbOrRef}; cutoff = 10) =
    Arblib.rel_accuracy_bits(x) < precision(x) - cutoff
iswide(x::Union{ArbSeries,AcbSeries}; cutoff = 10) = iswide(Arblib.ref(x, 0))
iswide(::Number; cutoff = 10) = false

"""
    mince(x::Arb, n::Integer)

Return a vector with `n` balls covering the ball `x`.
"""
function mince(x::Arb, n::Integer)
    balls = Vector{Arb}(undef, n)
    xₗ, xᵤ = Arblib.getinterval(Arb, x)
    dx = (xᵤ - xₗ) / n
    for i in eachindex(balls)
        yₗ = xₗ + (i - 1) * dx
        yᵤ = xₗ + i * dx
        balls[i] = Arb((yₗ, yᵤ))
    end

    return balls
end

# Conversion between Arb and Interval
Base.convert(::Type{Interval{T}}, x::Arb) where {T} =
    Interval{T}(getinterval(BigFloat, x)...)

Arblib.Arb(x::Interval{Float64}) =
    Arb((IntervalArithmetic.inf(x), IntervalArithmetic.sup(x)))


function arb_dot!(
    res::Arblib.ArbLike,
    ::Ptr{Nothing},
    subtract::Integer,
    x::Union{Arblib.ArbVectorLike,Ptr{Arblib.arb_struct}},
    xstep::Integer,
    y::Union{Arblib.ArbVectorLike,Ptr{Arblib.arb_struct}},
    ystep::Integer,
    len::Integer;
    prec::Integer = Arblib._precision(res),
)
    ccall(
        Arblib.@libarb(arb_dot),
        Nothing,
        (
            Ref{Arblib.arb_struct},
            Ptr{Nothing},
            Cint,
            Ptr{Arblib.arb_struct},
            Int,
            Ptr{Arblib.arb_struct},
            Int,
            Int,
            Int,
        ),
        res,
        C_NULL,
        subtract,
        x,
        xstep,
        y,
        ystep,
        len,
        prec,
    )
end


function arb_dot!(
    res::Arblib.ArbLike,
    s::Arblib.ArbLike,
    subtract::Integer,
    x::Union{Arblib.ArbVectorLike,Ptr{Arblib.arb_struct}},
    xstep::Integer,
    y::Union{Arblib.ArbVectorLike,Ptr{Arblib.arb_struct}},
    ystep::Integer,
    len::Integer;
    prec::Integer = Arblib._precision(res),
)
    ccall(
        Arblib.@libarb(arb_dot),
        Nothing,
        (
            Ref{Arblib.arb_struct},
            Ref{Arblib.arb_struct},
            Cint,
            Ptr{Arblib.arb_struct},
            Int,
            Ptr{Arblib.arb_struct},
            Int,
            Int,
            Int,
        ),
        res,
        s,
        subtract,
        x,
        xstep,
        y,
        ystep,
        len,
        prec,
    )
end
