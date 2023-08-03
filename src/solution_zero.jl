function _solve_capd(u0, (p, κ), ξ₀, ξ₁)
    a, b, α, β = u0

    d, ω, σ, ϵ, δ = p.d, p.ω, p.σ, p.ϵ, p.δ

    # FIXME: Handle this better
    IntervalArithmetic.setformat(:standard, sigfigs = 17)

    input_u0 = join(u0, "\n")
    input_params = join(Any[d, ω, σ, ϵ, δ, κ], "\n")
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
    _solve_zero_step(κ::Arb, μ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb; degree = 5)

Compute the solution at `ξ = ξ₁` using the Taylor expansion at `ξ =
0`.

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
function _solve_zero_step(κ::Arb, μ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb; degree = 20)
    # Compute expansion
    a, b = gl_taylor_expansion_real(
        SVector{2,NTuple{2,Arb}}((μ, 0), (0, 0)),
        zero(ξ₁),
        (p, κ);
        degree,
    )

    if p.σ == 1
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
                (1 + p.ϵ) * (
                    κ / (N + p.d) +
                    p.ω / ((N + 2) * (N + p.d)) +
                    (1 + p.δ) * (1 + 6M * C^3 / (N + p.d))
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

function _solve_zero_step(κ::T, μ::T, p::AbstractGLParams{T}, ξ₁::T; degree = 20) where {T}
    res = _solve_zero_step(
        convert(Arb, κ),
        convert(Arb, μ),
        gl_params(Arb, p),
        convert(Arb, ξ₁);
        degree,
    )

    return convert(SVector{4,T}, res)
end

"""
    solution_zero_capd(κ::T, μ::T, p::AbstractGLParams{T}, ξ₁::T) where {T}
    solution_zero_capd(κ::T, μ::T, p::AbstractGLParams{T}, ξ₀::T, ξ₁::T) where {T}

Integrate the system from `0` to `ξ₁` using the rigorous CAPD
integrator.

If `p.d != 1` the removable singularity at zero is handled using a
Taylor expansion at zero.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` and CAPD on `[ξ₀, ξ₁]`.
"""
function solution_zero_capd(κ::T, μ::T, p::AbstractGLParams{T}, ξ₁::T) where {T}
    if isone(p.d)
        u0 = SVector{4,Interval{Float64}}(μ, 0, 0, 0)
        ξ₀ = Interval{Float64}(0)
    else
        ξ₀ = convert(T, 1e-2)
        @assert ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        u0 = convert(SVector{4,Interval{Float64}}, _solve_zero_step(κ, μ, p, ξ₀))
        ξ₀ = Interval{Float64}(ξ₀)
    end

    p = gl_params(Interval{Float64}, p)
    κ = Interval{Float64}(κ)
    ξ₁ = Interval{Float64}(ξ₁)

    res = _solve_capd(u0, (p, κ), ξ₀, ξ₁)

    if T == Float64
        return IntervalArithmetic.mid.(res)
    else
        return convert.(T, res)
    end
end

function solution_zero_capd(κ::T, μ::T, p::AbstractGLParams{T}, ξ₀::T, ξ₁::T) where {T}
    @assert ξ₀ < ξ₁

    # Integrate system on [0, ξ₀] using Taylor expansion at zero
    u0 = convert(SVector{4,Interval{Float64}}, _solve_zero_step(κ, μ, p, ξ₀))

    # Integrate system on [ξ₀, ξ₁] using capd
    p = gl_params(Interval{Float64}, p)
    κ = Interval{Float64}(κ)
    ξ₀ = Interval{Float64}(ξ₀)
    ξ₁ = Interval{Float64}(ξ₁)

    res = _solve_capd(u0, (p, κ), ξ₀, ξ₁)

    if T == Float64
        return IntervalArithmetic.mid.(res)
    else
        return convert.(T, res)
    end
end

function solution_zero_float(
    κ::Float64,
    μ::Float64,
    p::AbstractGLParams{Float64},
    ξ₁::Float64,
)
    prob = ODEProblem(gl_equation_real, SVector(μ, 0, 0, 0), (0.0, ξ₁), (p, κ))

    sol = solve(prob, abstol = 1e-9, reltol = 1e-9)

    return sol[end]
end

function solution_zero_float(κ::Arb, μ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb)
    p = gl_params(Float64, p)
    ξ₁ = Float64(ξ₁)

    if iswide(κ)
        κs = collect(Float64.(getinterval(κ)))
    else
        κs = Float64(κ)
    end
    if iswide(μ)
        μs = collect(Float64.(getinterval(μ)))
    else
        μs = Float64(μ)
    end

    us = map(Iterators.product(κs, μs)) do (κ, μ)
        solution_zero_float(κ, μ, p, ξ₁)
    end

    return [Arb(extrema(getindex.(us, i))) for i in eachindex(us[begin])]
end
