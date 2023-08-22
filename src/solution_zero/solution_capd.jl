"""
    _solve_capd(
    u0::SVector{4,Interval{Float64}},
    κ::Interval{Float64},
    ξ₀::Interval{Float64},
    ξ₁::Interval{Float64},
    λ::AbstractGLParams{Interval{Float64}},
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
"""
function _solve_zero_capd(
    u0::SVector{4,Interval{Float64}},
    κ::Interval{Float64},
    ξ₀::Interval{Float64},
    ξ₁::Interval{Float64},
    λ::AbstractGLParams{Interval{Float64}},
)
    # FIXME: Handle this better
    IntervalArithmetic.setformat(:standard, sigfigs = 17)

    input_u0 = join(u0, "\n")
    input_params = join(Any[λ.d, λ.ω, λ.σ, λ.ϵ, λ.δ, κ], "\n")
    input_ξspan = join([ξ₀, ξ₁], "\n")

    input = join([input_u0, input_params, input_ξspan], "\n")

    # IMPROVE: Write directly to stdout of cmd instead of using echo
    program = pkgdir(@__MODULE__, "capd", "build", "ginzburg")
    cmd = pipeline(`echo $input`, `$program`)

    # FIXME: Avoid having to add this
    ENV["LD_LIBRARY_PATH"] = "/home/joeldahne/Programs/capd/lib/"

    output = readchomp(cmd)

    u1 = parse.(Interval{Float64}, split(output, "\n"))

    return u1
end

"""
TODO
"""
function _solve_zero_capd_jacobian(
    u0::SVector{4,Interval{Float64}},
    κ::Interval{Float64},
    ξ₀::Interval{Float64},
    ξ₁::Interval{Float64},
    λ::AbstractGLParams{Interval{Float64}},
)
    # FIXME: Handle this better
    IntervalArithmetic.setformat(:standard, sigfigs = 17)

    input_u0 = join(u0, "\n")
    input_params = join(Any[λ.d, λ.ω, λ.σ, λ.ϵ, λ.δ, κ], "\n")
    input_ξspan = join([ξ₀, ξ₁], "\n")

    input = join([input_u0, input_params, input_ξspan], "\n")

    # IMPROVE: Write directly to stdout of cmd instead of using echo
    program = pkgdir(@__MODULE__, "capd", "build", "ginzburg-variational")
    cmd = pipeline(`echo $input`, `$program`)

    # FIXME: Avoid having to add this
    ENV["LD_LIBRARY_PATH"] = "/home/joeldahne/Programs/capd/lib/"

    output = readchomp(cmd)

    res = parse.(Interval{Float64}, split(output, "\n"))

    u1 = res[1:4]
    du1_du01 = res[5:8]
    du1_dκ = res[9:12]

    return u1, du1_du01, du1_dκ
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
    if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        u0 = convert(SVector{4,Interval{Float64}}, _solve_zero_step(μ, κ, ξ₀, λ))
    else
        u0 = SVector{4,Interval{Float64}}(μ, 0, 0, 0)
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    κ = convert(Interval{Float64}, κ)
    ξ₀ = convert(Interval{Float64}, ξ₀)
    ξ₁ = convert(Interval{Float64}, ξ₁)
    λ = gl_params(Interval{Float64}, λ)

    res = _solve_zero_capd(u0, κ, ξ₀, ξ₁, λ)

    if T == Float64
        return IntervalArithmetic.mid.(res)
    else
        return convert.(T, res)
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
    if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        u0 = convert(SVector{4,Interval{Float64}}, _solve_zero_step(μ, κ, ξ₀, λ))
    else
        u0 = SVector{4,Interval{Float64}}(μ, 0, 0, 0)
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    κ = convert(Interval{Float64}, κ)
    ξ₀ = convert(Interval{Float64}, ξ₀)
    ξ₁ = convert(Interval{Float64}, ξ₁)
    λ = gl_params(Interval{Float64}, λ)

    res = _solve_zero_capd(u0, κ, ξ₀, ξ₁, λ)

    # TODO: Implement this!
    jacobian = SMatrix{4,2,Interval{Float64}}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    if T == Float64
        return IntervalArithmetic.mid.(res), IntervalArithmetic.mid.(jacobian)
    else
        return convert.(T, res), convert.(T, jacobian)
    end
end
