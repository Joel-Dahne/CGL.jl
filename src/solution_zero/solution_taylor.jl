"""
    _solution_zero_taylor_remainder(
        a::ArbSeries,
        b::ArbSeries,
        κ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term in
[`solution_zero_taylor`](@ref).

Let `N` be the degree of the arguments `a` and `b`. The remainder term
is bounded by finding `r` such that `abs(a[n])` and `abs(b[n])` are
bounded by `r^n` for `n > N`. The remainder term is then bounded as
```
abs(sum(a[n] * ξ₁^n for n = degree+1:Inf)) <= (r * ξ₁)^(degree + 1) / (1 - r * ξ₁)
```
and similarly for `b`. For the derivatives we instead get the bound
```
abs(sum(n * a[n] * ξ₁^(n - 1) for n = degree+1:Inf)) <=
    (r * ξ₁)^degree * (degree + 1 - degree * r * ξ₁) / (1 - r * ξ₁)^2
```
and the same for `b`.

For details on how we find `r` see lemma:tail-bound in the paper
(commit 095eee9).
"""
function _solution_zero_taylor_remainder(
    a::ArbSeries,
    b::ArbSeries,
    κ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    @assert Arblib.degree(a) == Arblib.degree(b)

    indeterminate_result = (indeterminate(κ), indeterminate(κ))

    isfinite(a) && isfinite(b) || return indeterminate_result

    N = Arblib.degree(a)

    if λ.σ == 1
        # r such that abs(a[n]), abs(b[n]) < r^n for n > N
        r = let
            # Value of r is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r = inv(16ξ₁)

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
                (1 + λ.ϵ) / (1 + λ.ϵ^2) * (
                    κ / (N + λ.d) +
                    λ.ω / ((N + 2) * (N + λ.d)) +
                    (1 + λ.δ) * (1 + 6M * C^3 / (N + λ.d))
                )

            D <= r^2 || return indeterminate_result

            r
        end

        @assert 0 < r * ξ₁ < 1

        remainder_bound = (r * ξ₁)^(N + 1) / (1 - r * ξ₁)
        remainder_derivative_bound = (r * ξ₁)^N * (N + 1 - N * r * ξ₁) / (1 - r * ξ₁)^2

        remainder = add_error(Arb(0), remainder_bound)
        remainder_derivative = add_error(Arb(0), remainder_derivative_bound)
    else
        @error "No implementation of remainder for σ != 1"

        remainder, remainder_derivative = zero(ξ₁), zero(ξ₁)
    end

    return remainder, remainder_derivative
end

"""
    _solution_zero_taylor_remainder_dμ(
        a::ArbSeries,
        b::ArbSeries,
        a_dμ::ArbSeries,
        b_dμ::ArbSeries,
        κ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term for the derivative w.r.t. μ
in [`solution_zero_jacobian_taylor`](@ref).

It works in the same way as [`_solution_zero_taylor_remainder`](@ref).

For details on how we find `r` see lemma:tail-bound-dmu in the paper.
"""
function _solution_zero_taylor_remainder_dμ(
    a::ArbSeries,
    b::ArbSeries,
    a_dμ::ArbSeries,
    b_dμ::ArbSeries,
    κ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    @assert Arblib.degree(a) ==
            Arblib.degree(b) ==
            Arblib.degree(a_dμ) ==
            Arblib.degree(b_dμ)

    indeterminate_result = (indeterminate(κ), indeterminate(κ))

    isfinite(a) && isfinite(b) && isfinite(a_dμ) && isfinite(b_dμ) ||
        return indeterminate_result

    N = Arblib.degree(a)

    if λ.σ == 1
        # r_μ such that abs(a_dμ[n]), abs(b_dμ[n]) < r_μ^n for n > N
        r_μ = let
            # Value of r_μ is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r_μ = inv(16ξ₁)

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
                (1 + λ.ϵ) / (1 + λ.ϵ^2) * (
                    κ / (N + λ.d) +
                    λ.ω / ((N + 2) * (N + λ.d)) +
                    3(1 + λ.δ) * (1 + 6M * C^3 / (N + λ.d))
                )

            D <= r_μ^2 || return indeterminate_result

            r_μ
        end

        @assert 0 < r_μ * ξ₁ < 1

        remainder_bound = (r_μ * ξ₁)^(N + 1) / (1 - r_μ * ξ₁)
        remainder_derivative_bound =
            (r_μ * ξ₁)^N * (N + 1 - N * r_μ * ξ₁) / (1 - r_μ * ξ₁)^2

        remainder = add_error(Arb(0), remainder_bound)
        remainder_derivative = add_error(Arb(0), remainder_derivative_bound)
    else
        @error "No implementation of remainder for σ != 1"

        remainder, remainder_derivative = zero(ξ₁), zero(ξ₁)
    end

    return remainder, remainder_derivative
end

"""
    _solution_zero_taylor_remainder_dκ(
        a::ArbSeries,
        b::ArbSeries,
        a_dκ::ArbSeries,
        b_dκ::ArbSeries,
        κ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term for the derivative w.r.t. κ
in [`solution_zero_jacobian_taylor`](@ref).

It works in the same way as [`_solution_zero_taylor_remainder`](@ref).

For details on how we find `r` see lemma:tail-bound-dkappa in the
paper.
"""
function _solution_zero_taylor_remainder_dκ(
    a::ArbSeries,
    b::ArbSeries,
    a_dκ::ArbSeries,
    b_dκ::ArbSeries,
    κ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    @assert Arblib.degree(a) ==
            Arblib.degree(b) ==
            Arblib.degree(a_dκ) ==
            Arblib.degree(b_dκ)

    indeterminate_result = (indeterminate(κ), indeterminate(κ))

    isfinite(a) && isfinite(b) && isfinite(a_dκ) && isfinite(b_dκ) ||
        return indeterminate_result

    N = Arblib.degree(a)

    if λ.σ == 1
        # r_κ such that abs(a_dκ[n]), abs(b_dκ[n]) < r_κ^n for n > N
        r_κ = let
            # Value of r_κ is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r_κ = inv(16ξ₁)

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
                (1 + λ.ϵ) / (1 + λ.ϵ^2) * (
                    (κ + 1) / (N + λ.d) +
                    λ.ω / ((N + 2) * (N + λ.d)) +
                    3(1 + λ.δ) * (1 + 6M * C^3 / (N + λ.d))
                )

            D <= r_κ^2 || return indeterminate_result

            r_κ
        end

        @assert 0 < r_κ * ξ₁ < 1

        remainder_bound = (r_κ * ξ₁)^(N + 1) / (1 - r_κ * ξ₁)
        remainder_derivative_bound =
            (r_κ * ξ₁)^N * (N + 1 - N * r_κ * ξ₁) / (1 - r_κ * ξ₁)^2

        remainder = add_error(Arb(0), remainder_bound)
        remainder_derivative = add_error(Arb(0), remainder_derivative_bound)
    else
        @error "No implementation of remainder for σ != 1"

        remainder, remainder_derivative = zero(ξ₁), zero(ξ₁)
    end

    return remainder, remainder_derivative
end

"""
    solution_zero_taylor(μ::T, κ::T, ξ₁::T, λ::CGLParams{T}; degree = 20) where {T}

Let `u = [a, b]` be a solution to [`ivp_zero_real`](@ref) This
function computes `[a(ξ₁), b(ξ₁), d(a)(ξ₁), d(b)(ξ₁)]` using the
Taylor expansion at `ξ = 0`.

This only works well for small values of `ξ₁` and is intended to be
used for handling the removable singularity at `ξ = 0`.
"""
function solution_zero_taylor(μ::Arb, κ::Arb, ξ₁::Arb, λ::CGLParams{Arb}; degree = 20)
    # Compute expansion
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        κ,
        zero(ξ₁),
        λ;
        degree,
    )

    remainder, remainder_derivative = _solution_zero_taylor_remainder(a, b, κ, ξ₁, λ)

    a0, a1 = Arblib.evaluate2(a, ξ₁)
    b0, b1 = Arblib.evaluate2(b, ξ₁)

    a0 += remainder
    a1 += remainder_derivative
    b0 += remainder
    b1 += remainder_derivative

    return SVector(a0, b0, a1, b1)
end

function solution_zero_taylor(μ::T, κ::T, ξ₁::T, λ::CGLParams{T}; degree = 20) where {T}
    u = solution_zero_taylor(
        convert(Arb, μ),
        convert(Arb, κ),
        convert(Arb, ξ₁),
        CGLParams{Arb}(λ);
        degree,
    )

    return convert(SVector{4,T}, u)
end

"""
    solution_zero_jacobian_taylor(μ::T, κ::T, ξ₁::T, λ::CGLParams{T}; degree = 20) where {T}

Let `u = [a, b]` be a solution to [`ivp_zero_real`](@ref) This
function computes `[a(ξ₁), b(ξ₁), d(a)(ξ₁), d(b)(ξ₁)]` using the
Taylor expansion at `ξ = 0`. It also computes the Jacobian w.r.t. `μ`
and `κ.

This only works well for small values of `ξ₁` and is intended to be
used for handling the removable singularity at `ξ = 0`.
"""
function solution_zero_jacobian_taylor(
    μ::Arb,
    κ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    degree = 20,
)
    # Compute expansion
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        κ,
        zero(ξ₁),
        λ;
        degree,
    )

    # Compute expansion of derivative w.r.t. μ
    a_dμ, b_dμ = cgl_equation_real_dμ_taylor(
        SVector{2,NTuple{2,Arb}}((1, 0), (0, 0)),
        a,
        b,
        κ,
        zero(ξ₁),
        λ;
        degree,
    )

    # Compute expansion of derivative w.r.t. κ
    a_dκ, b_dκ = cgl_equation_real_dκ_taylor(
        SVector{2,NTuple{2,Arb}}((0, 0), (0, 0)),
        a,
        b,
        κ,
        zero(ξ₁),
        λ;
        degree,
    )

    remainder, remainder_derivative = _solution_zero_taylor_remainder(a, b, κ, ξ₁, λ)
    remainder_dμ, remainder_derivative_dμ =
        _solution_zero_taylor_remainder_dμ(a, b, a_dμ, b_dμ, κ, ξ₁, λ)
    remainder_dκ, remainder_derivative_dκ =
        _solution_zero_taylor_remainder_dκ(a, b, a_dκ, b_dκ, κ, ξ₁, λ)

    a0, a1 = Arblib.evaluate2(a, ξ₁)
    b0, b1 = Arblib.evaluate2(b, ξ₁)
    a0_dμ, a1_dμ = Arblib.evaluate2(a_dμ, ξ₁)
    b0_dμ, b1_dμ = Arblib.evaluate2(b_dμ, ξ₁)
    a0_dκ, a1_dκ = Arblib.evaluate2(a_dκ, ξ₁)
    b0_dκ, b1_dκ = Arblib.evaluate2(b_dκ, ξ₁)

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

    u = SVector(a0, b0, a1, b1)

    J = SMatrix{4,2,Arb}(a0_dμ, b0_dμ, a1_dμ, b1_dμ, a0_dκ, b0_dκ, a1_dκ, b1_dκ)

    return u, J
end

function solution_zero_jacobian_taylor(
    μ::T,
    κ::T,
    ξ₁::T,
    λ::CGLParams{T};
    degree = 20,
) where {T}
    u, jacobian = solution_zero_jacobian_taylor(
        convert(Arb, μ),
        convert(Arb, κ),
        convert(Arb, ξ₁),
        CGLParams{Arb}(λ);
        degree,
    )

    return convert(SVector{4,T}, u), convert(SMatrix{4,2,T}, jacobian)
end
