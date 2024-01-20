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

# Conversion between Arb and BareInterval
Base.convert(::Type{BareInterval{T}}, x::Arb) where {T} =
    if isnan(x)
        # No nai constructor for BareInterval
        IntervalArithmetic.nai(T).bareinterval
    else
        bareinterval(T, getinterval(BigFloat, x)...)
    end

Base.convert(::Type{Arb}, x::BareInterval{Float64}) = Arb(x)

Arblib.Arb(x::BareInterval{Float64}) =
    if IntervalArithmetic.isempty_interval(x)
        indeterminate(Arb)
    else
        Arb((IntervalArithmetic.inf(x), IntervalArithmetic.sup(x)))
    end

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
        Arblib.@libflint(arb_dot),
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
        Arblib.@libflint(arb_dot),
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

    if sgn == 0 && !isinteger(y)
        # How many derivatives are well defined depends on the value
        # of y.

        # We need at most two derivatives in this case and therefore
        # only implement those. There could be more derivatives that
        # are finite but we don't need that.
        res[0] = abspow(x[0], y)
        if Arblib.degree(res) >= 1
            res[1] = y * abspow(x[0], y - 1) * x[1]
        end
        if Arblib.degree(res) >= 2
            res[2] =
                y * ((y - 1) * abspow(x[0], y - 2) * x[1]^2 + abspow(x[0], y - 1) * 2x[2]) /
                2
        end
        for i = 3:Arblib.degree(res)
            res[i] = indeterminate(Arb)
        end
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

# This is currently missing from Arblib due to parsing issues, we
# define it here instead.
Arblib.ArbCall.arbcall"void acb_rising2_ui(acb_t u, acb_t v, const acb_t x, ulong n, slong prec)"
