# This is currently missing from Arblib due to parsing issues, we
# define it here instead.
Arblib.ArbCall.arbcall"void acb_rising2_ui(acb_t u, acb_t v, const acb_t x, ulong n, slong prec)"

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
indeterminate(::Type{T}) where {T<:AbstractFloat} = convert(T, NaN)
indeterminate(::Type{Complex{T}}) where {T<:AbstractFloat} =
    convert(Complex{T}, complex(NaN, NaN))
indeterminate(x) = indeterminate(typeof(x))

"""
    iswide(x; cutoff = 10)

Return true if `x` is wide, in the meaning that the effective relative
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
abspow(x::Union{Arb,ArbSeries}, y::Arb) = abspow!(zero(x), x, y)

"""
    format_interval_precise(x; min_digits = 2)

Return a nicely formatted string representing the ball `x`.

Exact balls are printed exactly, non-exact balls are printed on the
form `1.234₅⁶`.

TODO: Implement support for ... in decimal expansion.
"""
function format_interval_precise(x::Arb; min_digits::Integer = 2)
    min_digits >= 1 || throw(ArgumentError("min_digits should be positive"))

    get_exact_string(y) =
        let digits = Arblib.digits_prec(precision(x)), res = Arblib.string(x; digits)
            @assert Arblib.isexact(y)
            while startswith(res, "[")
                digits *= 2
                res = string(y; digits)
                digits > 2^20 && error("$res") # Failsafe on large input
            end
            return res
        end

    if !isfinite(x)
        if isnan(x)
            res = "NaN"
        elseif Arblib.is_pos_inf(x)
            res = "\\infty"
        elseif Arblib.is_neg_inf(x)
            res = "-\\infty"
        else
            res = "[-\\infty, \\infty]"
        end
    elseif Arblib.isexact(x)
        res = get_exact_string(x)
    elseif Arblib.isnegative(x)
        res = "-" * format_interval_precise(-x; min_digits)
    elseif Arblib.ispositive(x)
        r = radius(Arf, x)
        m = Arblib.midref(x)
        ARF_PREC_EXACT = typemax(Int) # Flint constant

        upp = zero(x)
        exact = iszero(Arblib.add!(Arblib.midref(upp), m, r, prec = ARF_PREC_EXACT))
        @assert exact

        low = zero(x)
        exact = iszero(Arblib.sub!(Arblib.midref(low), m, r, prec = ARF_PREC_EXACT))
        @assert exact

        upp_string = get_exact_string(upp)
        low_string = get_exact_string(low)

        i = findfirst(1:min(length(low_string), length(upp_string))) do i
            low_string[i] != upp_string[i]
        end

        @assert !isnothing(i) # No common digits

        res_main = upp_string[1:i-1]

        if !contains(res_main, ".")
            # Currently only implement error after decimal point
            # TODO: Implement this
            @warn "Error before decimap point $x."
            return replace(string(x), "+/-" => "\\pm", r"e(.[0-9]*)" => s"\\cdot 10^{\1}")
        end

        if Arblib.ispositive(x)
            res_upp_start = upp_string[i:min(i + min_digits - 2, end)]
            res_low_start = low_string[i:min(i + min_digits - 2, end)]

            upp_string_remaining = upp_string[i+min_digits-1:end]
            low_string_remaining = low_string[i+min_digits-1:end]

            j = findfirst(!=(9), upp_string_remaining)

            if isnothing(j)
                res_upp_end = upp_string_remaining[1:end]
                res_low_end = low_string_remaining[1:end]
            else
                res_upp_end =
                    upp_string_remaining[1:j-1] *
                    string(parse(Int, upp_string_remaining[j]) + 1)
                res_low_end = low_string_remaining[1:min(j, end)]
            end

            res_superscript = res_upp_start * res_upp_end
            res_subscript = res_low_start * res_low_end
        elseif Arblib.negative(x)
            res_upp_start = upp_string[i:min(i + min_digits - 2, end)]
            res_low_start = low_string[i:min(i + min_digits - 2, end)]

            upp_string_remaining = upp_string[i+min_digits-1:end]
            low_string_remaining = low_string[i+min_digits-1:end]

            j = findfirst(!=(9), upp_string_remaining)

            if isnothing(j)
                res_upp_end = upp_string_remaining[1:end]
                res_low_end = low_string_remaining[1:end]
            else
                res_upp_end =
                    upp_string_remaining[1:j-1] *
                    string(parse(Int, upp_string_remaining[j]) + 1)
                res_low_end = low_string_remaining[1:min(j, end)]
            end

            res_superscript = res_upp_start * res_upp_end
            res_subscript = res_low_start * res_low_end
        end

        res = res_main * "_{" * res_subscript * "}^{" * res_superscript * "}"
    else
        # TODO: Implement this
        @warn "Contains zero $x."
        return replace(string(x), "+/-" => "\\pm", r"e(.[0-9]*)" => s"\\cdot 10^{\1}")
    end

    return res
end

function format_interval_precise(z::Acb; min_digits::Integer = 2)
    re_str = format_interval_precise(real(z); min_digits)
    im_str = format_interval_precise(imag(z); min_digits)

    if startswith(im_str, "-")
        return re_str * " - " * im_str[2:end] * " i"
    else
        return re_str * " + " * im_str * " i"
    end
end
