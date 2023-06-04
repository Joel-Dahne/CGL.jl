function _solve_capd(u0, (p, κ), ξ₀, ξ₁)
    a, b, α, β = u0

    d, ω, σ, ϵ, δ = p.d, p.ω, p.σ, p.ϵ, p.δ

    # FIXME: Handle this better
    IntervalArithmetic.displaymode(sigfigs = 17)

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

function solution_zero_capd(κ::T, μ::T, p::AbstractGLParams{T}, ξ₁::T) where {T}
    u0 = SVector{4,Interval{Float64}}(μ, 0, 0, 0)
    p = gl_params(Interval{Float64}, p)
    κ = Interval{Float64}(κ)
    ξ₀ = Interval{Float64}(0)
    ξ₁ = Interval{Float64}(ξ₁)

    res = _solve_capd(u0, (p, κ), ξ₀, ξ₁)

    if T == Float64
        return IntervalArithmetic.mid.(res)
    else
        return convert.(T, res)
    end
end

function solution_zero_float(κ::T, μ::T, p::AbstractGLParams{T}, ξ₁::T) where {T}
    prob = ODEProblem(
        gl_equation_real,
        SVector{4,Float64}(μ, 0, 0, 0),
        (0.0, Float64(ξ₁)),
        (gl_params(Float64, p), Float64(κ)),
    )

    sol = solve(prob, abstol = 1e-9, reltol = 1e-9)

    return sol[end]
end
