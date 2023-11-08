"""
    _solve_zero_capd(
        u0::SVector{4,Interval{Float64}},
        κ::Interval{Float64},
        ξ₀::Interval{Float64},
        ξ₁::Interval{Float64},
        λ::AbstractGLParams{Interval{Float64}},
        output_jacobian::Union{Val{false},Val{true}} = Val{false}(),
    )

Let `u = [a, b, α, β]` be a solution to
[`ivp_zero_real_system`](@ref), but with initial values
```
a(ξ₀) = u0[1]
b(ξ₀) = u0[2]
α(ξ₀) = u0[3]
β(ξ₀) = u0[4]
```
This function computes `u(ξ₁)`.

If `output_jacobian = Val{true}()` it also computes the Jacobian
w.r.t. `[a, b, α, β, κ]`.

- **IMPROVE:** For thin values of `κ` and `λ.ϵ` the performance could
  be improved when computing `u`. For thin values of `λ.ϵ` it could be
  improved also when computing the Jacobian.
"""
function _solve_zero_capd(
    u0::SVector{4,Interval{Float64}},
    κ::Interval{Float64},
    ξ₀::Interval{Float64},
    ξ₁::Interval{Float64},
    λ::AbstractGLParams{Interval{Float64}},
    output_jacobian::Union{Val{false},Val{true}} = Val{false}(),
)
    input_u0 = ""
    for x in u0
        input_u0 *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_params = "$(λ.d)\n"
    for x in [λ.ω, λ.σ, λ.ϵ, λ.δ, κ]
        input_params *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_ξspan = ""
    for x in [ξ₀, ξ₁]
        input_ξspan *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_output_jacobian = ifelse(output_jacobian isa Val{false}, "0", "1")

    input = join([input_u0, input_params, input_ξspan, input_output_jacobian])

    # IMPROVE: Write directly to stdout of cmd instead of using echo
    program = pkgdir(@__MODULE__, "capd", "build", "ginzburg")
    cmd = pipeline(`echo $input`, `$program`)

    # FIXME: Avoid having to add this
    ENV["LD_LIBRARY_PATH"] = "/home/joeldahne/Programs/capd/lib/"

    output = readchomp(cmd)

    res = parse.(Interval{Float64}, split(output, "\n"))::Vector{Interval{Float64}}

    u1 = SVector{4,Interval{Float64}}(res[1:4])

    if output_jacobian isa Val{false}
        return u1
    else
        J = SMatrix{4,5,Interval{Float64}}(res[5:end])
        return u1, J
    end
end

"""
    _solve_zero_step(μ::T, κ::T, ξ₁::T, λ::AbstractGLParams{T}; degree = 20) where {T}

Let `u = [a, b]` be a solution to [`ivp_zero_real`](@ref) This
function computes `[a(ξ₁), b(ξ₁), d(a)(ξ₁), d(b)(ξ₁)]` using the
Taylor expansion at `ξ = 0`.

This only works well for small values of `ξ₁` and is intended to be
used for handling the removable singularity at `ξ = 0`.

# Bounding remainder term
The remainder term is bounded by finding `r` such that `abs(a[n])` and
`abs(b[n])` are bounded by `r^n` for `n > degree`. The remainder term
for is then bounded as
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
function _solve_zero_step(μ::Arb, κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb}; degree = 20)
    # Compute expansion
    a, b = gl_equation_real_taylor_expansion(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        κ,
        zero(ξ₁),
        λ;
        degree,
    )

    if λ.σ == 1
        # r such that abs(a[n]), abs(b[n]) < r^n for n > degree
        r = let N = degree
            # Value of r is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r = inv(4ξ₁)

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
            @assert all(n -> abs(a[n]) <= C * r^n, 0:M-1)
            @assert all(n -> abs(a[n]) <= r^n, M:degree)
            @assert all(n -> abs(b[n]) <= C * r^n, 0:M-1)
            @assert all(n -> abs(b[n]) <= r^n, M:degree)

            # This is needed for the lemma to apply
            3M < N || error("3M < N not satisfied, M = $M, N = $N")

            D =
                (1 + λ.ϵ) * (
                    κ / (N + λ.d) +
                    λ.ω / ((N + 2) * (N + λ.d)) +
                    (1 + λ.δ) * (1 + 6M * C^3 / (N + λ.d))
                )

            D <= r^2 || error("D < r^2 not satisfied, r = $r, D = $D")

            r
        end

        @assert 0 < r * ξ₁ < 1

        remainder, remainder_derivative = let N = degree
            remainder_bound = (r * ξ₁)^(N + 1) / (1 - r * ξ₁)
            remainder_derivative_bound =
                (r * ξ₁)^N * (N + 1 - N * r * ξ₁) / (1 - r * ξ₁)^2

            add_error(Arb(0), remainder_bound),
            add_error(Arb(0), remainder_derivative_bound)
        end
    else
        @error "No implementation of remainder for σ != 1"

        remainder, remainder_derivative = zero(ξ₁), zero(ξ₁)
    end

    a0, a1 = Arblib.evaluate2(a, ξ₁)
    b0, b1 = Arblib.evaluate2(b, ξ₁)

    a0 += remainder
    a1 += remainder_derivative
    b0 += remainder
    b1 += remainder_derivative

    return SVector(a0, b0, a1, b1)
end

function _solve_zero_step(μ::T, κ::T, ξ₁::T, λ::AbstractGLParams{T}; degree = 20) where {T}
    res = _solve_zero_step(
        convert(Arb, μ),
        convert(Arb, κ),
        convert(Arb, ξ₁),
        gl_params(Arb, λ);
        degree,
    )

    return convert(SVector{4,T}, res)
end

"""
    _solve_zero_jacobian_step(μ::T, κ::T, ξ₁::T, λ::AbstractGLParams{T}; degree = 20) where {T}

Let `u = [a, b]` be a solution to [`ivp_zero_real`](@ref) This
function computes `[a(ξ₁), b(ξ₁), d(a)(ξ₁), d(b)(ξ₁)]` using the
Taylor expansion at `ξ = 0`. It also computes the Jacobian w.r.t. `μ`
and `κ.

This only works well for small values of `ξ₁` and is intended to be
used for handling the removable singularity at `ξ = 0`.
"""
function _solve_zero_jacobian_step(
    μ::Arb,
    κ::Arb,
    ξ₁::Arb,
    λ::AbstractGLParams{Arb};
    degree = 20,
)
    # Compute expansion
    a, b = gl_equation_real_taylor_expansion(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        κ,
        zero(ξ₁),
        λ;
        degree,
    )

    if λ.σ == 1
        # r such that abs(a[n]), abs(b[n]) < r^n for n > degree
        r = let N = degree
            # Value of r is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r = inv(4ξ₁)

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
            @assert all(n -> abs(a[n]) <= C * r^n, 0:M-1)
            @assert all(n -> abs(a[n]) <= r^n, M:degree)
            @assert all(n -> abs(b[n]) <= C * r^n, 0:M-1)
            @assert all(n -> abs(b[n]) <= r^n, M:degree)

            # This is needed for the lemma to apply
            3M < N || error("3M < N not satisfied, M = $M, N = $N")

            D =
                (1 + λ.ϵ) * (
                    κ / (N + λ.d) +
                    λ.ω / ((N + 2) * (N + λ.d)) +
                    (1 + λ.δ) * (1 + 6M * C^3 / (N + λ.d))
                )

            D <= r^2 || error("D < r^2 not satisfied, r = $r, D = $D")

            r
        end

        @assert 0 < r * ξ₁ < 1

        remainder, remainder_derivative = let N = degree
            remainder_bound = (r * ξ₁)^(N + 1) / (1 - r * ξ₁)
            remainder_derivative_bound =
                (r * ξ₁)^N * (N + 1 - N * r * ξ₁) / (1 - r * ξ₁)^2

            add_error(Arb(0), remainder_bound),
            add_error(Arb(0), remainder_derivative_bound)
        end
    else
        @error "No implementation of remainder for σ != 1"

        remainder, remainder_derivative = zero(ξ₁), zero(ξ₁)
    end

    a0, a1 = Arblib.evaluate2(a, ξ₁)
    b0, b1 = Arblib.evaluate2(b, ξ₁)

    a0 += remainder
    a1 += remainder_derivative
    b0 += remainder
    b1 += remainder_derivative

    u1 = SVector(a0, b0, a1, b1)

    # TODO: For now we approximate the Jacobian using the Float64
    # version
    #jacobian = SMatrix{4,2,Arb}(1, 0, 0, 0, 0, 0, 0, 0)
    jacobian = convert(
        SMatrix{4,2,Arb},
        solution_zero_jacobian_float(
            Float64(μ),
            Float64(κ),
            Float64(ξ₁),
            gl_params(Float64, λ),
        )[2],
    )

    return u1, jacobian
end

function _solve_zero_jacobian_step(
    μ::T,
    κ::T,
    ξ₁::T,
    λ::AbstractGLParams{T};
    degree = 20,
) where {T}
    u1, jacobian = _solve_zero_jacobian_step(
        convert(Arb, μ),
        convert(Arb, κ),
        convert(Arb, ξ₁),
        gl_params(Arb, λ);
        degree,
    )

    return convert(SVector{4,T}, u1), convert(SMatrix{4,2,T}, jacobian)
end

"""
    solution_zero_capd(μ::T, κ::T, ξ₁::T, λ::AbstractGLParams{T}) where {T}
    solution_zero_capd(μ::T, κ::T, ξ₀::T, ξ₁::T, λ::AbstractGLParams{T}) where {T}

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes `u(ξ₁)`.

The solution is computed using the rigorous CAPD integrator.

If `p.d != 1` the removable singularity at zero is handled using a
Taylor expansion at zero.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` and CAPD on `[ξ₀, ξ₁]`.
"""
function solution_zero_capd(μ::T, κ::T, ξ₁::T, λ::AbstractGLParams{T}) where {T}
    if isone(λ.d)
        ξ₀ = zero(ξ₁)
    else
        ξ₀ = convert(T, 1e-2)
    end

    return solution_zero_capd(μ, κ, ξ₀, ξ₁, λ)
end

function solution_zero_capd(μ::T, κ::T, ξ₀::T, ξ₁::T, λ::AbstractGLParams{T}) where {T}
    S = Interval{Float64}

    u0 = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        convert(SVector{4,S}, _solve_zero_step(μ, κ, ξ₀, λ))
    else
        SVector{4,S}(interval(Float64, μ), interval(0.0), interval(0.0), interval(0.0))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    u1 = let
        κ = convert(S, κ)
        ξ₀ = convert(S, ξ₀)
        ξ₁ = convert(S, ξ₁)
        λ = gl_params(S, λ)

        _solve_zero_capd(u0, κ, ξ₀, ξ₁, λ)
    end

    if T == Float64
        return IntervalArithmetic.mid.(u1)
    else
        return convert.(T, u1)
    end
end

"""
    solution_zero_jacobian_capd(μ::T, κ::T, ξ₁::T, λ::AbstractGLParams{T}) where {T}
    solution_zero_jacobian_capd(μ::T, κ::T, ξ₀::T, ξ₁::T, λ::AbstractGLParams{T}) where {T}

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes `u(ξ₁)` as well as the Jacobian w.r.t. `μ` and
`κ`.

The solution is computed using the rigorous CAPD integrator.

If `p.d != 1` the removable singularity at zero is handled using a
Taylor expansion at zero.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` and CAPD on `[ξ₀, ξ₁]`.
"""
function solution_zero_jacobian_capd(μ::T, κ::T, ξ₁::T, λ::AbstractGLParams{T}) where {T}
    if isone(λ.d)
        ξ₀ = zero(ξ₁)
    else
        ξ₀ = convert(T, 1e-2)
    end

    return solution_zero_jacobian_capd(μ, κ, ξ₀, ξ₁, λ)
end

function solution_zero_jacobian_capd(
    μ::T,
    κ::T,
    ξ₀::T,
    ξ₁::T,
    λ::AbstractGLParams{T},
) where {T}
    S = Interval{Float64}

    u0, J1 = let
        if !iszero(ξ₀)
            @assert 0 < ξ₀ < ξ₁
            # Integrate system on [0, ξ₀] using Taylor expansion at zero
            u0, J1 = _solve_zero_jacobian_step(μ, κ, ξ₀, λ)
            u0 = convert(SVector{4,S}, u0)
            J1 = convert(SMatrix{4,2,S}, J1)
        else
            u0 = SVector{4,S}(interval(Float64, μ), interval(0.0), interval(0.0), interval(0.0))
            # Empty integration so the only non-zero derivative is the
            # one of u0[1] w.r.t. μ, which is 1.
            J1 = SMatrix{4,2,S}(
                interval(1.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
            )
        end
        # J1 now contains derivatives of [u0[1], u0[2], u0[3], u0[4]].
        # We want to add a row [0, 1] for the derivative of κ.
        u0, vcat(J1, SMatrix{1,2,S}(interval(0.0), interval(1.0)))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    u1, J2 = let
        κ = convert(S, κ)
        ξ₀ = convert(S, ξ₀)
        ξ₁ = convert(S, ξ₁)
        λ = gl_params(S, λ)

        _solve_zero_capd(u0, κ, ξ₀, ξ₁, λ, Val{true}())
    end

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    J = J2 * J1

    if T == Float64
        return IntervalArithmetic.mid.(u1), IntervalArithmetic.mid.(J)
    else
        return convert.(T, u1), convert.(T, J)
    end
end
