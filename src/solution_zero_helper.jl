"""
    _u1nu2n!(u1n::Arb, u2n::Arb, a2b2::ArbRefVector, a2b2σ::ArbRefVector, a::ArbSeries, b::ArbSeries, σ::Arb, n::Integer)

This function is used for computing the `n`th Taylor coefficient of
the two functions
```
u1 = (a^2 + b^2)^σ * a
u2 = (a^2 + b^2)^σ * b
```

The arguments `a` and `b` should hold all Taylor coefficients up to
the `n`th one. The arguments `a2b2σ` and `a2b2` should hold the Taylor
coefficients of `(a^2 + b^2)` and `(a^2 + b^2)^σ` up to, but not
including, the `n`th one, they should have space for the `n`th one.
"""
function _u1nu2n!(
    u1n::Arb,
    u2n::Arb,
    a2b2::ArbRefVector,
    a2b2σ::ArbRefVector,
    a::ArbSeries,
    b::ArbSeries,
    σ::Arb,
    n::Integer,
)
    # Pointer to first and last entries in a and b to use
    a_first = a.poly.arb_poly.coeffs
    b_first = b.poly.arb_poly.coeffs
    a_last = a.poly.arb_poly.coeffs + n * sizeof(Arblib.arb_struct)
    b_last = b.poly.arb_poly.coeffs + n * sizeof(Arblib.arb_struct)

    # Set coefficient for a2b2 vector
    arb_dot!(a2b2[n+1], C_NULL, 0, a_first, 1, a_last, -1, n + 1)
    arb_dot!(a2b2[n+1], a2b2[n+1], 0, b_first, 1, b_last, -1, n + 1)

    # IMPROVE: Split this into its separate parts and optimize them
    # individually
    Arblib.pow_arb_series!(a2b2σ, a2b2, n + 1, σ, n + 1)

    arb_dot!(u1n, C_NULL, 0, a2b2σ, 1, a_last, -1, n + 1)
    arb_dot!(u2n, C_NULL, 0, a2b2σ, 1, b_last, -1, n + 1)

    return nothing
end
