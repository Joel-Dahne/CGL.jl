"""
    _u1nu2n!(u1n::Arb, u2n::Arb, a2b2::ArbRefVector, a2b2σ::ArbRefVector, a::ArbSeries, b::ArbSeries, σ::Arb, n::Integer)

This function is used for computing the `n`th Taylor coefficient of
the two functions
```
u1 = (a^2 + b^2)^σ * a
u2 = (a^2 + b^2)^σ * b
```

The computation of `(a^2 + b^2)^σ` is done as
```
(a^2 + b^2)^σ = exp(σ * log(a^2 + b^2))
```
For the logarithm it uses that
```
log(x) = Arblib.integral(Arblib.derivative(x) / x)
```
with the constant of integration set to `log(x[0])`.

For the computations it uses several intermediate terms. These are
reused between computations and therefore taken in as arguments. It
uses the following arguments.
- `a` and `b`: Should hold all Taylor coefficients up to the `n`th
  one. These are not mutated.
- `a2b2, σ_log_a2b2, a2b2σ`: Correspond to `a^2 + b^2`, `σ * log(a^2 +
  b^2)` and `(a^2 + b^2)^σ` respectively. Should hold the Taylor
  coefficients up to, but not including, the `n`th one. Should have
  space for the `n`th coefficient, which is set during the
  calculations.
- `a2b2_diff, a2b2_inv, a2b2_div`: Correspond to
  `Arblib.derivative(a^2 + b^2)`, `inv(a^2 + b^2)` and
  `Arblib.derivative(a^2 + b^2) / (a^2 + b^2)` respectively. Should
  hold all the Taylor coefficients up to, but not including, the `n -
  1`th one. Should have space for the `n - 1`th coefficient, which is
  set during the calculations.
"""
function _u1nu2n!(
    u1n::Arb,
    u2n::Arb,
    a2b2::ArbRefVector,
    a2b2_diff::ArbRefVector,
    a2b2_inv::ArbRefVector,
    a2b2_div::ArbRefVector,
    σ_log_a2b2::ArbRefVector,
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

    # Compute coefficient for  a^2 + b^2
    arb_dot!(a2b2[n+1], C_NULL, 0, a_first, 1, a_last, -1, n + 1)
    arb_dot!(a2b2[n+1], a2b2[n+1], 0, b_first, 1, b_last, -1, n + 1)

    # Compute coefficient for σ * log(a^2 + b^2)
    if iszero(n)
        Arblib.log!(σ_log_a2b2[1], a2b2[1])
        Arblib.mul!(σ_log_a2b2[1], σ_log_a2b2[1], σ)
    else
        # Compute coefficient of Arblib.derivative(a^2 + b^2)
        Arblib.mul!(a2b2_diff[n], a2b2[n+1], n)

        # Compute coefficient of inv(a^2 + b^2)
        # IMPROVE: This currently computes the full series every time
        Arblib.inv_series!(a2b2_inv, a2b2, n + 1, n)

        # Compute coefficient of Arblib.derivative(a^2 + b^2) * inv(a^2 + b^2)
        a2b2_inv_last = a2b2_inv.arb_vec.entries + (n - 1) * sizeof(Arblib.arb_struct)
        arb_dot!(a2b2_div[n], C_NULL, 0, a2b2_diff, 1, a2b2_inv_last, -1, n)

        # Compute coefficient of
        # Arblib.integral(Arblib.derivative(a^2 + b^2) / (a^2 + b^2))
        Arblib.div!(σ_log_a2b2[n+1], a2b2_div[n], n)

        Arblib.mul!(σ_log_a2b2[n+1], σ_log_a2b2[n+1], σ)
    end

    # Compute coefficient for exp(σ * log(a^2 + b^2))
    # IMPROVE: This currently computes the full series every time
    Arblib.exp_series!(a2b2σ, σ_log_a2b2, n + 1, n + 1)

    arb_dot!(u1n, C_NULL, 0, a2b2σ, 1, a_last, -1, n + 1)
    arb_dot!(u2n, C_NULL, 0, a2b2σ, 1, b_last, -1, n + 1)

    return nothing
end
