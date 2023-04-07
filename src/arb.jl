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
