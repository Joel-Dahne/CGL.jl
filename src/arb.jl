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
    interval(T, getinterval(BigFloat, x)...)

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

"""
    abspow!(res, x, y)

Inplace version of [`abspow`](@ref).
"""
function abspow!(res::Arb, x::Arblib.ArbOrRef, y::Arb)
    iszero(y) && return Arblib.one!(res)

    if iszero(x)
        Arblib.contains_negative(y) && return Arblib.indeterminate!(res)
        Arblib.ispositive(y) && return Arblib.zero!(res)
        return Arblib.unit_interval!(res)
    end

    if Arblib.contains_zero(x)
        Arblib.contains_negative(y) && return Arblib.indeterminate!(res)
        upper = abs_ubound(Arb, x) # One extra allocation
        Arblib.pow!(upper, upper, y)
        Arblib.zero!(res)
        return Arblib.union!(res, res, upper)
    end

    if res === y
        # In this case we need an extra allocation to not overwrite y
        y = copy(y)
    end
    Arblib.abs!(res, x)
    return Arblib.pow!(res, res, y)
end

function abspow!(res::ArbSeries, x::ArbSeries, y::Arb)
    Arblib.degree(res) == Arblib.degree(x) ||
        throw(ArgumentError("res and x should have the same degree"))

    sgn = Arblib.sgn_nonzero(Arblib.ref(x, 0))

    if sgn == 0
        # We don't have to be that careful with allocations here.

        # All non-constant terms are indeterminate, the constant term
        # is given by abs(x[0])^y
        res[0] = abspow(x[0], y)
        for i = 1:Arblib.degree(res)
            res[i] = indeterminate(Arb)
        end

        return res
    elseif sgn < 0
        Arblib.neg!(res, x)
        Arblib.pow_arb_series!(res, res, y, length(res))
    else
        Arblib.pow_arb_series!(res, x, y, length(res))
    end

    return res
end

"""
    abspow(x, y)

Compute `abs(x)^y `in a way that works if `x` overlaps with zero.
"""
abspow(x::Arb, y::Arb) = abspow!(zero(x), x, y)
abspow(x::ArbSeries, y::Arb) = abspow!(zero(x), x, y)
