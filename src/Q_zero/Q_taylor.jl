"""
    _Q_zero_taylor_remainder(
        a::ArbSeries,
        b::ArbSeries,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        λ::CGLParams{Arb},
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

For details on how we find `r` see lemma:tail-bound in the paper
(commit 095eee9).
"""
function _Q_zero_taylor_remainder(
    a::ArbSeries,
    b::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb},
)
    @assert Arblib.degree(a) == Arblib.degree(b)

    indeterminate_result = (indeterminate(κ), indeterminate(κ), indeterminate(κ))

    isfinite(a) && isfinite(b) || return indeterminate_result

    (; d, ω, σ, δ) = λ

    N = Arblib.degree(a)

    if σ == 1
        # r such that abs(a[n]), abs(b[n]) < r^n for n > N
        r = let
            # Value of r is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r = inv(16ξ₀)

            # Find C such that abs(a[n]), abs(b[n]) < C * r^n for 0 <= k <= N
            C =
                1.01max(
                    maximum(n -> abs(a[n] / r^n), 0:N),
                    maximum(n -> abs(b[n] / r^n), 0:N),
                )

            # Find M such that abs(a[n]), abs(b[n]) <= r^n for M <= n <= N
            M = let M = findlast(n -> !(abs(a[n]) <= r^n && abs(b[n]) <= r^n), 0:N)
                isnothing(M) ? 0 : M
            end

            # Verify that r, C1 and M satisfy the requirements
            all(n -> abs(a[n]) <= C * r^n, 0:M-1) || return indeterminate_result
            all(n -> abs(a[n]) <= r^n, M:N) || return indeterminate_result
            all(n -> abs(b[n]) <= C * r^n, 0:M-1) || return indeterminate_result
            all(n -> abs(b[n]) <= r^n, M:N) || return indeterminate_result

            # This is needed for the lemma to apply
            3M < N || return indeterminate_result

            D =
                (1 + abs(ϵ)) / (1 + ϵ^2) * (
                    abs(κ) / (N + d) +
                    abs(ω) / ((N + 2) * (N + d)) +
                    (1 + abs(δ)) * (1 + 6M * C^3 / (N + d))
                )

            D <= r^2 || return indeterminate_result

            r
        end

        @assert 0 < r * ξ₀ < 1

        remainder_bound = (r * ξ₀)^(N + 1) / (1 - r * ξ₀)
        remainder_derivative_bound = (r * ξ₀)^N * (N + 1 - N * r * ξ₀) / (1 - r * ξ₀)^2
        remainder_derivative2_bound =
            (r * ξ₀)^(N - 1) * (N + N^2 + (2 - 2N^2) * r * ξ₀ - (N - N^2) * (r * ξ₀)^2) /
            (1 - r * ξ₀)^3

        remainder = add_error(Arb(0), remainder_bound)
        remainder_derivative = add_error(Arb(0), remainder_derivative_bound)
        remainder_derivative2 = add_error(Arb(0), remainder_derivative2_bound)
    else
        error("No implementation of remainder for σ != 1")
    end

    return remainder, remainder_derivative, remainder_derivative2
end

"""
    _Q_zero_taylor_remainder_dμ(
        a::ArbSeries,
        b::ArbSeries,
        a_dμ::ArbSeries,
        b_dμ::ArbSeries,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term for the derivative w.r.t. μ
in [`Q_zero_jacobian_taylor`](@ref).

It works in the same way as [`_Q_zero_taylor_remainder`](@ref), except
it doesn't return a bound for the second derivative.

For details on how we find `r` see lemma:tail-bound-dmu in the paper.
"""
function _Q_zero_taylor_remainder_dμ(
    a::ArbSeries,
    b::ArbSeries,
    a_dμ::ArbSeries,
    b_dμ::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb},
)
    @assert Arblib.degree(a) ==
            Arblib.degree(b) ==
            Arblib.degree(a_dμ) ==
            Arblib.degree(b_dμ)

    indeterminate_result = (indeterminate(κ), indeterminate(κ))

    isfinite(a) && isfinite(b) && isfinite(a_dμ) && isfinite(b_dμ) ||
        return indeterminate_result

    (; d, ω, σ, δ) = λ

    N = Arblib.degree(a)

    if σ == 1
        # r_μ such that abs(a_dμ[n]), abs(b_dμ[n]) < r_μ^n for n > N
        r_μ = let
            # Value of r_μ is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r_μ = inv(16ξ₀)

            # Find C such that
            # abs(a[n]), abs(b[n]), abs(a_dμ[n]), abs(b_dμ[n]) < C * r^n
            # for 0 <= k <= N
            C =
                1.01max(
                    maximum(n -> abs(a[n] / r_μ^n), 0:N),
                    maximum(n -> abs(b[n] / r_μ^n), 0:N),
                    maximum(n -> abs(a_dμ[n] / r_μ^n), 0:N),
                    maximum(n -> abs(b_dμ[n] / r_μ^n), 0:N),
                )

            # Find M such that
            # abs(a[n]), abs(b[n]), abs(a_dμ[n]), abs(b_dμ[n]) < r^n
            # for M <= n <= N
            M =
                let M = findlast(
                        n -> !(
                            abs(a[n]) <= r_μ^n &&
                            abs(b[n]) <= r_μ^n &&
                            abs(a_dμ[n]) <= r_μ^n &&
                            abs(b_dμ[n]) <= r_μ^n
                        ),
                        0:N,
                    )
                    isnothing(M) ? 0 : M
                end

            # Verify that r, C1 and M satisfy the requirements
            all(n -> abs(a[n]) <= C * r_μ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(a[n]) <= r_μ^n, M:N) || return indeterminate_result
            all(n -> abs(b[n]) <= C * r_μ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(b[n]) <= r_μ^n, M:N) || return indeterminate_result
            all(n -> abs(a_dμ[n]) <= C * r_μ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(a_dμ[n]) <= r_μ^n, M:N) || return indeterminate_result
            all(n -> abs(b_dμ[n]) <= C * r_μ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(b_dμ[n]) <= r_μ^n, M:N) || return indeterminate_result

            # This is needed for the lemma to apply
            3M < N || return indeterminate_result

            D =
                (1 + abs(ϵ)) / (1 + ϵ^2) * (
                    abs(κ) / (N + d) +
                    abs(ω) / ((N + 2) * (N + d)) +
                    3(1 + abs(δ)) * (1 + 6M * C^3 / (N + d))
                )

            D <= r_μ^2 || return indeterminate_result

            r_μ
        end

        @assert 0 < r_μ * ξ₀ < 1

        remainder_bound = (r_μ * ξ₀)^(N + 1) / (1 - r_μ * ξ₀)
        remainder_derivative_bound =
            (r_μ * ξ₀)^N * (N + 1 - N * r_μ * ξ₀) / (1 - r_μ * ξ₀)^2

        remainder = add_error(Arb(0), remainder_bound)
        remainder_derivative = add_error(Arb(0), remainder_derivative_bound)
    else
        @error "No implementation of remainder for σ != 1"

        remainder, remainder_derivative = zero(ξ₀), zero(ξ₀)
    end

    return remainder, remainder_derivative
end

"""
    _Q_zero_taylor_remainder_dκ(
        a::ArbSeries,
        b::ArbSeries,
        a_dκ::ArbSeries,
        b_dκ::ArbSeries,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term for the derivative w.r.t. κ
in [`Q_zero_jacobian_taylor`](@ref).

It works in the same way as [`_Q_zero_taylor_remainder`](@ref), except
it doesn't return a bound for the second derivative.

For details on how we find `r` see lemma:tail-bound-dkappa in the
paper.
"""
function _Q_zero_taylor_remainder_dκ(
    a::ArbSeries,
    b::ArbSeries,
    a_dκ::ArbSeries,
    b_dκ::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb},
)
    @assert Arblib.degree(a) ==
            Arblib.degree(b) ==
            Arblib.degree(a_dκ) ==
            Arblib.degree(b_dκ)

    indeterminate_result = (indeterminate(κ), indeterminate(κ))

    isfinite(a) && isfinite(b) && isfinite(a_dκ) && isfinite(b_dκ) ||
        return indeterminate_result

    (; d, ω, σ, δ) = λ

    N = Arblib.degree(a)

    if σ == 1
        # r_κ such that abs(a_dκ[n]), abs(b_dκ[n]) < r_κ^n for n > N
        r_κ = let
            # Value of r_κ is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r_κ = inv(16ξ₀)

            # Find C such that
            # abs(a[n]), abs(b[n]), abs(a_dκ[n]), abs(b_dκ[n]) < C * r^n
            # for 0 <= k <= N
            C =
                1.01max(
                    maximum(n -> abs(a[n] / r_κ^n), 0:N),
                    maximum(n -> abs(b[n] / r_κ^n), 0:N),
                    maximum(n -> abs(a_dκ[n] / r_κ^n), 0:N),
                    maximum(n -> abs(b_dκ[n] / r_κ^n), 0:N),
                )

            # Find M such that
            # abs(a[n]), abs(b[n]), abs(a_dκ[n]), abs(b_dκ[n]) < r^n
            # for M <= n <= N
            M =
                let M = findlast(
                        n -> !(
                            abs(a[n]) <= r_κ^n &&
                            abs(b[n]) <= r_κ^n &&
                            abs(a_dκ[n]) <= r_κ^n &&
                            abs(b_dκ[n]) <= r_κ^n
                        ),
                        0:N,
                    )
                    isnothing(M) ? 0 : M
                end

            # Verify that r, C1 and M satisfy the requirements
            all(n -> abs(a[n]) <= C * r_κ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(a[n]) <= r_κ^n, M:N) || return indeterminate_result
            all(n -> abs(b[n]) <= C * r_κ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(b[n]) <= r_κ^n, M:N) || return indeterminate_result
            all(n -> abs(a_dκ[n]) <= C * r_κ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(a_dκ[n]) <= r_κ^n, M:N) || return indeterminate_result
            all(n -> abs(b_dκ[n]) <= C * r_κ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(b_dκ[n]) <= r_κ^n, M:N) || return indeterminate_result

            # This is needed for the lemma to apply
            3M < N || return indeterminate_result

            D =
                (1 + abs(ϵ)) / (1 + ϵ^2) * (
                    (abs(κ) + 1) / (N + d) +
                    abs(ω) / ((N + 2) * (N + d)) +
                    3(1 + abs(δ)) * (1 + 6M * C^3 / (N + d))
                )

            D <= r_κ^2 || return indeterminate_result

            r_κ
        end

        @assert 0 < r_κ * ξ₀ < 1

        remainder_bound = (r_κ * ξ₀)^(N + 1) / (1 - r_κ * ξ₀)
        remainder_derivative_bound =
            (r_κ * ξ₀)^N * (N + 1 - N * r_κ * ξ₀) / (1 - r_κ * ξ₀)^2

        remainder = add_error(Arb(0), remainder_bound)
        remainder_derivative = add_error(Arb(0), remainder_derivative_bound)
    else
        @error "No implementation of remainder for σ != 1"

        remainder, remainder_derivative = zero(ξ₀), zero(ξ₀)
    end

    return remainder, remainder_derivative
end

"""
    _Q_zero_taylor_remainder_dϵ(
        a::ArbSeries,
        b::ArbSeries,
        a_dμ::ArbSeries,
        b_dμ::ArbSeries,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term for the derivative w.r.t. ϵ
in [`Q_zero_jacobian_taylor`](@ref).

It works in the same way as [`_Q_zero_taylor_remainder`](@ref), except
it doesn't return a bound for the second derivative.

For details on how we find `r` see lemma:tail-bound-depsilon in the
paper.
"""
function _Q_zero_taylor_remainder_dϵ(
    a::ArbSeries,
    b::ArbSeries,
    a_dϵ::ArbSeries,
    b_dϵ::ArbSeries,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb},
)
    @assert Arblib.degree(a) ==
            Arblib.degree(b) ==
            Arblib.degree(a_dϵ) ==
            Arblib.degree(b_dϵ)

    indeterminate_result = (indeterminate(κ), indeterminate(κ))

    isfinite(a) && isfinite(b) && isfinite(a_dϵ) && isfinite(b_dϵ) ||
        return indeterminate_result

    (; d, ω, σ, δ) = λ

    N = Arblib.degree(a)

    if σ == 1
        # r_ϵ such that abs(a_dϵ[n]), abs(b_dϵ[n]) < r_ϵ^n for n > N
        r_ϵ = let
            # Value of r_ϵ is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r_ϵ = inv(16ξ₀)

            # Find C such that
            # abs(a[n]), abs(b[n]), abs(a_dϵ[n]), abs(b_dϵ[n]) < C * r^n
            # for 0 <= k <= N
            C =
                1.01max(
                    maximum(n -> abs(a[n] / r_ϵ^n), 0:N),
                    maximum(n -> abs(b[n] / r_ϵ^n), 0:N),
                    maximum(n -> abs(a_dϵ[n] / r_ϵ^n), 0:N),
                    maximum(n -> abs(b_dϵ[n] / r_ϵ^n), 0:N),
                )

            # Find M such that
            # abs(a[n]), abs(b[n]), abs(a_dϵ[n]), abs(b_dϵ[n]) < r^n
            # for M <= n <= N
            M =
                let M = findlast(
                        n -> !(
                            abs(a[n]) <= r_ϵ^n &&
                            abs(b[n]) <= r_ϵ^n &&
                            abs(a_dϵ[n]) <= r_ϵ^n &&
                            abs(b_dϵ[n]) <= r_ϵ^n
                        ),
                        0:N,
                    )
                    isnothing(M) ? 0 : M
                end

            # Verify that r, C1 and M satisfy the requirements
            all(n -> abs(a[n]) <= C * r_ϵ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(a[n]) <= r_ϵ^n, M:N) || return indeterminate_result
            all(n -> abs(b[n]) <= C * r_ϵ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(b[n]) <= r_ϵ^n, M:N) || return indeterminate_result
            all(n -> abs(a_dϵ[n]) <= C * r_ϵ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(a_dϵ[n]) <= r_ϵ^n, M:N) || return indeterminate_result
            all(n -> abs(b_dϵ[n]) <= C * r_ϵ^n, 0:M-1) || return indeterminate_result
            all(n -> abs(b_dϵ[n]) <= r_ϵ^n, M:N) || return indeterminate_result

            # This is needed for the lemma to apply
            3M < N || return indeterminate_result

            D =
                (1 + abs(ϵ)) / (1 + ϵ^2) * (
                    1 +
                    abs(κ) / (N + d) +
                    abs(ω) / ((N + 2) * (N + d)) +
                    3(1 + abs(δ)) * (1 + 6M * C^3 / (N + d))
                )

            D <= r_ϵ^2 || return indeterminate_result

            r_ϵ
        end

        @assert 0 < r_ϵ * ξ₀ < 1

        remainder_bound = (r_ϵ * ξ₀)^(N + 1) / (1 - r_ϵ * ξ₀)
        remainder_derivative_bound =
            (r_ϵ * ξ₀)^N * (N + 1 - N * r_ϵ * ξ₀) / (1 - r_ϵ * ξ₀)^2

        remainder = add_error(Arb(0), remainder_bound)
        remainder_derivative = add_error(Arb(0), remainder_derivative_bound)
    else
        @error "No implementation of remainder for σ != 1"

        remainder, remainder_derivative = zero(ξ₀), zero(ξ₀)
    end

    return remainder, remainder_derivative
end

"""
    Q_zero_taylor(μ, κ, ϵ, ξ₀, λ::CGLParams; degree = 20, enclose_curve = Val{false}())

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
    λ::CGLParams{Arb};
    degree = 20,
    enclose_curve::Union{Val{false},Val{true}} = Val{false}(),
)
    # Compute expansion
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        κ,
        ϵ,
        zero(ξ₀),
        λ;
        degree,
    )

    remainder, remainder_derivative, remainder_derivative2 =
        _Q_zero_taylor_remainder(a, b, κ, ϵ, ξ₀, λ)

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

function Q_zero_taylor(μ::T, κ::T, ϵ::T, ξ₀::T, λ::CGLParams{T}; degree = 20) where {T}
    Q = Q_zero_taylor(Arb(μ), Arb(κ), Arb(ϵ), Arb(ξ₀), CGLParams{Arb}(λ); degree)

    return T.(Q)
end

"""
    Q_zero_jacobian_kappa_taylor(μ, κ, ϵ, ξ₀, λ::CGLParams; degree = 20)

This function computes the Jacobian of [`Q_zero_taylor`](@ref) w.r.t.
the parameters `μ` and `κ`. It also returns the result of
[`Q_zero_taylor`](@ref).
"""
function Q_zero_jacobian_kappa_taylor(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree = 20,
)
    # Compute expansion
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        κ,
        ϵ,
        zero(ξ₀),
        λ;
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
        λ;
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
        λ;
        degree,
    )

    remainder, remainder_derivative, _ = _Q_zero_taylor_remainder(a, b, κ, ϵ, ξ₀, λ)
    remainder_dμ, remainder_derivative_dμ =
        _Q_zero_taylor_remainder_dμ(a, b, a_dμ, b_dμ, κ, ϵ, ξ₀, λ)
    remainder_dκ, remainder_derivative_dκ =
        _Q_zero_taylor_remainder_dκ(a, b, a_dκ, b_dκ, κ, ϵ, ξ₀, λ)

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

function Q_zero_jacobian_kappa_taylor(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₀::T,
    λ::CGLParams{T};
    degree = 20,
) where {T}
    Q, J = Q_zero_jacobian_taylor(
        Arb.(μ),
        Arb.(κ),
        Arb.(ϵ),
        Arb.(ξ₀),
        CGLParams{Arb}(λ);
        degree,
    )

    return T.(Q), T.(J)
end

"""
    Q_zero_jacobian_epsilon_taylor(μ, κ, ϵ, ξ₀, λ::CGLParams; degree = 20)

This function computes the Jacobian of [`Q_zero_taylor`](@ref) w.r.t.
the parameters `μ` and `ϵ`. It also returns the result of
[`Q_zero_taylor`](@ref).
"""
function Q_zero_jacobian_epsilon_taylor(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree = 20,
)
    # Compute expansion
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        κ,
        ϵ,
        zero(ξ₀),
        λ;
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
        λ;
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
        λ;
        degree,
    )

    remainder, remainder_derivative, _ = _Q_zero_taylor_remainder(a, b, κ, ϵ, ξ₀, λ)
    remainder_dμ, remainder_derivative_dμ =
        _Q_zero_taylor_remainder_dμ(a, b, a_dμ, b_dμ, κ, ϵ, ξ₀, λ)
    remainder_dϵ, remainder_derivative_dϵ =
        _Q_zero_taylor_remainder_dϵ(a, b, a_dϵ, b_dϵ, κ, ϵ, ξ₀, λ)

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

function Q_zero_jacobian_epsilon_taylor(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₀::T,
    λ::CGLParams{T};
    degree = 20,
) where {T}
    Q, J = Q_zero_jacobian_taylor(
        Arb.(μ),
        Arb.(κ),
        Arb.(ϵ),
        Arb.(ξ₀),
        CGLParams{Arb}(λ);
        degree,
    )

    return T.(Q), T.(J)
end
