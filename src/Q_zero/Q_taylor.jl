"""
    _Q_zero_taylor_remainder_check_conditions(M, N, C, r, a, b, κ, ϵ, Λ)

Check that `M`,`N`, `C` and `r` satisfy the conditions of Lemma
REF(lemma:tail-bound).
"""
function _Q_zero_taylor_remainder_check_conditions(
    M::Int,
    N::Int,
    C::Arb,
    r::Arb,
    a::ArbSeries,
    b::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    Λ::CGLParams{Arb},
)
    # Ensure a and b have coefficient up to N computed
    @assert Arblib.degree(a) == Arblib.degree(b) == N

    (; d, ω, σ, δ) = Λ
    @assert isone(σ) # The lemma is only valid for σ = 1

    iseven(N) || return false
    3M < N || return false

    all(n -> abs(a[n]) <= C * r^n, 0:(M-1)) || return false
    all(n -> abs(b[n]) <= C * r^n, 0:(M-1)) || return false
    all(n -> abs(a[n]) <= r^n, M:N) || return false
    all(n -> abs(b[n]) <= r^n, M:N) || return false

    m = (M + 1) ÷ 2 # This is the same as ceil(M / 2)

    D =
        (1 + abs(ϵ)) / (1 + ϵ^2) * (
            abs(κ) / (N + d) +
            abs(ω) / ((N + 2) * (N + d)) +
            2(1 + abs(δ)) * (
                1 // 8 +
                1 // 2N +
                3m * C * ((1 + 3 // N) // 2(N + d)) +
                3m^2 * C^2 / ((N + 2) * (N + d))
            )
        )

    return D <= r^2
end

"""
    _Q_zero_taylor_remainder(
        a::ArbSeries,
        b::ArbSeries,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        Λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term in [`Q_zero_taylor`](@ref).

Let `N` be the degree of the arguments `a` and `b`. The remainder term
is bounded by finding `r` such that `abs(a[n])` and `abs(b[n])` are
bounded by `r^n` for `n > N`. The remainder term is then bounded as
```
abs(sum(a[n] * ξ₀^n for n = N+1:Inf)) <= (r * ξ₀)^(N + 1) / (1 - r * ξ₀)
```
and similarly for `b`. For the derivatives we instead get the bound
```
abs(sum(n * a[n] * ξ₀^(n - 1) for n = N+1:Inf)) <=
    (r * ξ₀)^N * (N + 1 - N * r * ξ₀) / (1 - r * ξ₀)^2
```
and the same for `b`. For the second derivative we get the bound
```
abs(sum(n * (n - 1) * a[n] * ξ₀^(n - 2) for n = N+1:Inf)) <=
    (r * ξ₀)^(N - 1) * (N + N^2 + (2 - 2N^2) * r * ξ₀ - (N - N^2) * (r * ξ₀)^2) / (1 - r * ξ₀)^3
```
The bound for the second derivative is used in
[`verification_monotonicity`](@ref).

Given `r`, the proof that `abs(a[n])` and `abs(b[n])` are bounded by
`r^N` for `n > N` is based on Lemma REF(lemma:tail-bound).
"""
function _Q_zero_taylor_remainder(
    a::ArbSeries,
    b::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb},
)
    @assert Arblib.degree(a) == Arblib.degree(b)

    (; d, ω, σ, δ) = Λ
    N = Arblib.degree(a)

    isone(σ) || error("No implementation of remainder for σ != 1")

    isfinite(a) && isfinite(b) ||
        return indeterminate(Arb), indeterminate(Arb), indeterminate(Arb)

    # Value of r is a tuning parameter. Lower value gives tighter
    # enclosures but makes it harder to verify the requirements.
    r = inv(16ξ₀)
    @assert 0 < r * ξ₀ < 1 # Sanity check

    # Find C such that abs(a[n]), abs(b[n]) < C * r^n for 0 <= k <= N
    C = ubound(
        Arb,
        1.01max(maximum(n -> abs(a[n] / r^n), 0:N), maximum(n -> abs(b[n] / r^n), 0:N)),
    )

    # Find M such that abs(a[n]), abs(b[n]) <= r^n for M <= n <= N
    M = let M = findlast(n -> !(abs(a[n]) <= r^n && abs(b[n]) <= r^n), 0:N)
        isnothing(M) ? 0 : M
    end

    # Check that the conditions of the lemma are satisfied
    if !_Q_zero_taylor_remainder_check_conditions(M, N, C, r, a, b, κ, ϵ, Λ)
        return indeterminate(Arb), indeterminate(Arb), indeterminate(Arb)
    end

    remainder_bound = (r * ξ₀)^(N + 1) / (1 - r * ξ₀)
    remainder_derivative_bound = (r * ξ₀)^N * (N + 1 - N * r * ξ₀) / (1 - r * ξ₀)^2
    remainder_derivative2_bound =
        (r * ξ₀)^(N - 1) * (N + N^2 + (2 - 2N^2) * r * ξ₀ - (N - N^2) * (r * ξ₀)^2) /
        (1 - r * ξ₀)^3

    remainder = add_error(Arb(0), remainder_bound)
    remainder_derivative = add_error(Arb(0), remainder_derivative_bound)
    remainder_derivative2 = add_error(Arb(0), remainder_derivative2_bound)

    return remainder, remainder_derivative, remainder_derivative2
end

"""
    _Q_zero_taylor_remainder_dμ(
        a_dμ::ArbSeries,
        b_dμ::ArbSeries,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        Λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term for the derivative w.r.t. μ
in [`Q_zero_jacobian_kappa_taylor`](@ref) and [`Q_zero_jacobian_epsilon_taylor`](@ref).

As discussed in the paper, we can use exactly the same lemma as for
`a` and `b` (without the derivatives). This is hence just a wrapper of
[`_Q_zero_taylor_remainder`](@ref). Contrary to
`_Q_zero_taylor_remainder`, it does however not return the bound for
the second derivative (we don't need it).
"""
_Q_zero_taylor_remainder_dμ(
    a_dμ::ArbSeries,
    b_dμ::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb},
) = _Q_zero_taylor_remainder(a_dμ, b_dμ, κ, ϵ, ξ₀, Λ)[1:2]

"""
    _Q_zero_taylor_remainder_dκ(
        a::ArbSeries,
        b::ArbSeries,
        a_dκ::ArbSeries,
        b_dκ::ArbSeries,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        Λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term for the derivative w.r.t. κ
in [`Q_zero_jacobian_kappa_taylor`](@ref).

It works in the same way as [`_Q_zero_taylor_remainder`](@ref), except
it doesn't return a bound for the second derivative.

See Lemma REF(lemma:tail-bound-dkappa) in the paper for details on how
the tail is bounded.
"""
function _Q_zero_taylor_remainder_dκ(
    a::ArbSeries,
    b::ArbSeries,
    a_dκ::ArbSeries,
    b_dκ::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb},
)
    @assert Arblib.degree(a) ==
            Arblib.degree(b) ==
            Arblib.degree(a_dκ) ==
            Arblib.degree(b_dκ)

    (; d, ω, σ, δ) = Λ
    N = Arblib.degree(a)

    isone(σ) || error("No implementation of remainder for σ != 1")

    isfinite(a) && isfinite(b) && isfinite(a_dκ) && isfinite(b_dκ) ||
        return indeterminate(Arb), indeterminate(Arb)

    # Value of r is a tuning parameter. Lower value gives tighter
    # enclosures but makes it harder to verify the requirements.
    r = inv(16ξ₀)
    @assert 0 < r * ξ₀ < 1 # Sanity check

    # Find C such that
    # abs(a[n]), abs(b[n]), abs(a_dκ[n]), abs(b_dκ[n]) < C * r^n
    # for 0 <= k <= N
    C = ubound(
        Arb,
        1.01max(
            maximum(n -> abs(a[n] / r^n), 0:N),
            maximum(n -> abs(b[n] / r^n), 0:N),
            maximum(n -> abs(a_dκ[n] / r^n), 0:N),
            maximum(n -> abs(b_dκ[n] / r^n), 0:N),
        ),
    )

    # Find M such that
    # abs(a[n]), abs(b[n]), abs(a_dκ[n]), abs(b_dκ[n]) < r^n
    # for M <= n <= N
    M =
        let M = findlast(
                n -> !(
                    abs(a[n]) <= r^n &&
                    abs(b[n]) <= r^n &&
                    abs(a_dκ[n]) <= r^n &&
                    abs(b_dκ[n]) <= r^n
                ),
                0:N,
            )
            isnothing(M) ? 0 : M
        end

    # Check that the conditions of the Lemma
    # REF(lemma:tail-bound) are satisfied (they are also a
    # requirement for Lemma REF(lemma:tail-bound-dkappa)).
    if !_Q_zero_taylor_remainder_check_conditions(M, N, C, r, a, b, κ, ϵ, Λ)
        return indeterminate(Arb), indeterminate(Arb)
    end

    # Check that the conditions for Lemma REF(lemma:tail-bound-dkappa) are satisfied
    ok = true
    ok &= all(n -> abs(a_dκ[n]) <= C * r^n, 0:(M-1))
    ok &= all(n -> abs(b_dκ[n]) <= C * r^n, 0:(M-1))
    ok &= all(n -> abs(a_dκ[n]) <= r^n, M:N)
    ok &= all(n -> abs(b_dκ[n]) <= r^n, M:N)

    ok || return indeterminate(Arb), indeterminate(Arb)

    m = (M + 1) ÷ 2 # This is the same as ceil(M / 2)

    D =
        (1 + abs(ϵ)) / (1 + ϵ^2) * (
            (abs(κ) + 1) / (N + d) +
            abs(ω) / ((N + 2) * (N + d)) +
            6(1 + abs(δ)) * (
                1 // 8 +
                1 // 2N +
                3m * C * ((1 + 3 // N) // 2(N + d)) +
                3m^2 * C^2 / ((N + 2) * (N + d))
            )
        )

    D <= r^2 || return indeterminate(Arb), indeterminate(Arb)

    remainder_bound = (r * ξ₀)^(N + 1) / (1 - r * ξ₀)
    remainder_derivative_bound = (r * ξ₀)^N * (N + 1 - N * r * ξ₀) / (1 - r * ξ₀)^2

    remainder = add_error(Arb(0), remainder_bound)
    remainder_derivative = add_error(Arb(0), remainder_derivative_bound)

    return remainder, remainder_derivative
end

"""
    _Q_zero_taylor_remainder_dϵ(
        a::ArbSeries,
        b::ArbSeries,
        a_dϵ::ArbSeries,
        b_dϵ::ArbSeries,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        Λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term for the derivative w.r.t. ϵ
in [`Q_zero_jacobian_epsilon_taylor`](@ref).

It works in the same way as [`_Q_zero_taylor_remainder`](@ref), except
it doesn't return a bound for the second derivative.

See Lemma REF(lemma:tail-bound-depsilon) in the paper for details on
how the tail is bounded.
"""
function _Q_zero_taylor_remainder_dϵ(
    a::ArbSeries,
    b::ArbSeries,
    a_dϵ::ArbSeries,
    b_dϵ::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb},
)
    @assert Arblib.degree(a) ==
            Arblib.degree(b) ==
            Arblib.degree(a_dϵ) ==
            Arblib.degree(b_dϵ)

    (; d, ω, σ, δ) = Λ
    N = Arblib.degree(a)

    isone(σ) || error("No implementation of remainder for σ != 1")

    isfinite(a) && isfinite(b) && isfinite(a_dϵ) && isfinite(b_dϵ) ||
        return indeterminate(Arb), indeterminate(Arb)

    # Value of r is a tuning parameter. Lower value gives tighter
    # enclosures but makes it harder to verify the requirements.
    r = inv(16ξ₀)
    @assert 0 < r * ξ₀ < 1 # Sanity check

    # Find C such that
    # abs(a[n]), abs(b[n]), abs(a_dϵ[n]), abs(b_dϵ[n]) < C * r^n
    # for 0 <= k <= N
    C = ubound(
        Arb,
        1.01max(
            maximum(n -> abs(a[n] / r^n), 0:N),
            maximum(n -> abs(b[n] / r^n), 0:N),
            maximum(n -> abs(a_dϵ[n] / r^n), 0:N),
            maximum(n -> abs(b_dϵ[n] / r^n), 0:N),
        ),
    )

    # Find M such that
    # abs(a[n]), abs(b[n]), abs(a_dϵ[n]), abs(b_dϵ[n]) < r^n
    # for M <= n <= N
    M =
        let M = findlast(
                n -> !(
                    abs(a[n]) <= r^n &&
                    abs(b[n]) <= r^n &&
                    abs(a_dϵ[n]) <= r^n &&
                    abs(b_dϵ[n]) <= r^n
                ),
                0:N,
            )
            isnothing(M) ? 0 : M
        end

    # Check that the conditions of the Lemma
    # REF(lemma:tail-bound) are satisfied (they are also a
    # requirement for Lemma REF(lemma:tail-bound-depsilon)).
    if !_Q_zero_taylor_remainder_check_conditions(M, N, C, r, a, b, κ, ϵ, Λ)
        return indeterminate(Arb), indeterminate(Arb)
    end

    # Check that the conditions for Lemma REF(lemma:tail-bound-depsilon) are satisfied
    ok = true
    ok &= all(n -> abs(a_dϵ[n]) <= C * r^n, 0:(M-1))
    ok &= all(n -> abs(b_dϵ[n]) <= C * r^n, 0:(M-1))
    ok &= all(n -> abs(a_dϵ[n]) <= r^n, M:N)
    ok &= all(n -> abs(b_dϵ[n]) <= r^n, M:N)

    ok || return indeterminate(Arb), indeterminate(Arb)

    m = (M + 1) ÷ 2 # This is the same as ceil(M / 2)

    D =
        (1 + abs(ϵ)) / (1 + ϵ^2) * (
            1 +
            abs(κ) / (N + d) +
            abs(ω) / ((N + 2) * (N + d)) +
            6(1 + abs(δ)) * (
                1 // 8 +
                1 // 2N +
                3m * C * ((1 + 3 // N) // 2(N + d)) +
                3m^2 * C^2 / ((N + 2) * (N + d))
            )
        )

    D <= r^2 || return indeterminate(Arb), indeterminate(Arb)

    remainder_bound = (r * ξ₀)^(N + 1) / (1 - r * ξ₀)
    remainder_derivative_bound = (r * ξ₀)^N * (N + 1 - N * r * ξ₀) / (1 - r * ξ₀)^2

    remainder = add_error(Arb(0), remainder_bound)
    remainder_derivative = add_error(Arb(0), remainder_derivative_bound)

    return remainder, remainder_derivative
end

"""
    Q_zero_taylor(μ, κ, ϵ, ξ₀, Λ::CGLParams; degree = 20, enclose_curve = Val{false}())

Compute the solution to the ODE on the interval ``[0, ξ₀]``. Returns a
vector with four real values, the first two are the real and imaginary
values at `ξ₀` and the second two are their derivatives.

The solution is computed using a Taylor expansion at `ξ = 0`. This
only works well for small values of `ξ₀` and is intended to be used
for handling the removable singularity at `ξ = 0`.

If `enclose_curve` is set to `Val{true}()` then don't return the value
at `ξ₀`, but an enclosure valid for `0 <= ξ <= ξ₀`. In this case it
also returns a second value containing enclosures of the real and
imaginary parts of the second order derivatives.
"""
function Q_zero_taylor(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb};
    degree = 20,
    enclose_curve::Union{Val{false},Val{true}} = Val{false}(),
)
    # Compute expansion
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        κ,
        ϵ,
        zero(ξ₀),
        Λ;
        degree,
    )

    remainder, remainder_derivative, remainder_derivative2 =
        _Q_zero_taylor_remainder(a, b, κ, ϵ, ξ₀, Λ)

    if enclose_curve isa Val{true}
        a0, a1 = Arblib.evaluate2(a, Arb((0, ξ₀)))
        b0, b1 = Arblib.evaluate2(b, Arb((0, ξ₀)))
    else
        a0, a1 = Arblib.evaluate2(a, ξ₀)
        b0, b1 = Arblib.evaluate2(b, ξ₀)
    end

    a0 += remainder
    a1 += remainder_derivative
    b0 += remainder
    b1 += remainder_derivative

    if enclose_curve isa Val{true}
        a2 = Arblib.derivative(a, 2)(Arb((0, ξ₀)))
        b2 = Arblib.derivative(b, 2)(Arb((0, ξ₀)))
        a2 += remainder_derivative2
        b2 += remainder_derivative2

        return SVector(a0, b0, a1, b1), SVector(a2, b2)
    else
        return SVector(a0, b0, a1, b1)
    end
end

"""
    Q_zero_jacobian_kappa_taylor(μ, κ, ϵ, ξ₀, Λ::CGLParams; degree = 20)

This function computes the Jacobian of [`Q_zero_taylor`](@ref) w.r.t.
the parameters `μ` and `κ`. It also returns the result of
[`Q_zero_taylor`](@ref).
"""
function Q_zero_jacobian_kappa_taylor(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb};
    degree = 20,
)
    # Compute expansion
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        κ,
        ϵ,
        zero(ξ₀),
        Λ;
        degree,
    )

    # Compute expansion of derivative w.r.t. μ
    a_dμ, b_dμ = cgl_equation_real_dμ_taylor(
        SVector{2,NTuple{2,Arb}}((1, 0), (0, 0)),
        a,
        b,
        κ,
        ϵ,
        zero(ξ₀),
        Λ;
        degree,
    )

    # Compute expansion of derivative w.r.t. κ
    a_dκ, b_dκ = cgl_equation_real_dκ_taylor(
        SVector{2,NTuple{2,Arb}}((0, 0), (0, 0)),
        a,
        b,
        κ,
        ϵ,
        zero(ξ₀),
        Λ;
        degree,
    )

    remainder, remainder_derivative, _ = _Q_zero_taylor_remainder(a, b, κ, ϵ, ξ₀, Λ)
    remainder_dμ, remainder_derivative_dμ =
        _Q_zero_taylor_remainder_dμ(a_dμ, b_dμ, κ, ϵ, ξ₀, Λ)
    remainder_dκ, remainder_derivative_dκ =
        _Q_zero_taylor_remainder_dκ(a, b, a_dκ, b_dκ, κ, ϵ, ξ₀, Λ)

    a0, a1 = Arblib.evaluate2(a, ξ₀)
    b0, b1 = Arblib.evaluate2(b, ξ₀)
    a0_dμ, a1_dμ = Arblib.evaluate2(a_dμ, ξ₀)
    b0_dμ, b1_dμ = Arblib.evaluate2(b_dμ, ξ₀)
    a0_dκ, a1_dκ = Arblib.evaluate2(a_dκ, ξ₀)
    b0_dκ, b1_dκ = Arblib.evaluate2(b_dκ, ξ₀)

    a0 += remainder
    a1 += remainder_derivative
    b0 += remainder
    b1 += remainder_derivative
    a0_dμ += remainder_dμ
    a1_dμ += remainder_derivative_dμ
    b0_dμ += remainder_dμ
    b1_dμ += remainder_derivative_dμ
    a0_dκ += remainder_dκ
    a1_dκ += remainder_derivative_dκ
    b0_dκ += remainder_dκ
    b1_dκ += remainder_derivative_dκ

    Q = SVector(a0, b0, a1, b1)

    J = SMatrix{4,2}(a0_dμ, b0_dμ, a1_dμ, b1_dμ, a0_dκ, b0_dκ, a1_dκ, b1_dκ)

    return Q, J
end

"""
    Q_zero_jacobian_epsilon_taylor(μ, κ, ϵ, ξ₀, Λ::CGLParams; degree = 20)

This function computes the Jacobian of [`Q_zero_taylor`](@ref) w.r.t.
the parameters `μ` and `ϵ`. It also returns the result of
[`Q_zero_taylor`](@ref).
"""
function Q_zero_jacobian_epsilon_taylor(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb};
    degree = 20,
)
    # Compute expansion
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        κ,
        ϵ,
        zero(ξ₀),
        Λ;
        degree,
    )

    # Compute expansion of derivative w.r.t. μ
    a_dμ, b_dμ = cgl_equation_real_dμ_taylor(
        SVector{2,NTuple{2,Arb}}((1, 0), (0, 0)),
        a,
        b,
        κ,
        ϵ,
        zero(ξ₀),
        Λ;
        degree,
    )

    # Compute expansion of derivative w.r.t. ϵ
    a_dϵ, b_dϵ = cgl_equation_real_dϵ_taylor(
        SVector{2,NTuple{2,Arb}}((0, 0), (0, 0)),
        a,
        b,
        κ,
        ϵ,
        zero(ξ₀),
        Λ;
        degree,
    )

    remainder, remainder_derivative, _ = _Q_zero_taylor_remainder(a, b, κ, ϵ, ξ₀, Λ)
    remainder_dμ, remainder_derivative_dμ =
        _Q_zero_taylor_remainder_dμ(a_dμ, b_dμ, κ, ϵ, ξ₀, Λ)
    remainder_dϵ, remainder_derivative_dϵ =
        _Q_zero_taylor_remainder_dϵ(a, b, a_dϵ, b_dϵ, κ, ϵ, ξ₀, Λ)

    a0, a1 = Arblib.evaluate2(a, ξ₀)
    b0, b1 = Arblib.evaluate2(b, ξ₀)
    a0_dμ, a1_dμ = Arblib.evaluate2(a_dμ, ξ₀)
    b0_dμ, b1_dμ = Arblib.evaluate2(b_dμ, ξ₀)
    a0_dϵ, a1_dϵ = Arblib.evaluate2(a_dϵ, ξ₀)
    b0_dϵ, b1_dϵ = Arblib.evaluate2(b_dϵ, ξ₀)

    a0 += remainder
    a1 += remainder_derivative
    b0 += remainder
    b1 += remainder_derivative
    a0_dμ += remainder_dμ
    a1_dμ += remainder_derivative_dμ
    b0_dμ += remainder_dμ
    b1_dμ += remainder_derivative_dμ
    a0_dϵ += remainder_dϵ
    a1_dϵ += remainder_derivative_dϵ
    b0_dϵ += remainder_dϵ
    b1_dϵ += remainder_derivative_dϵ

    Q = SVector(a0, b0, a1, b1)

    J = SMatrix{4,2}(a0_dμ, b0_dμ, a1_dμ, b1_dμ, a0_dϵ, b0_dϵ, a1_dϵ, b1_dϵ)

    return Q, J
end
