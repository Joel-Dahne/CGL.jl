using Dates

include("helper.jl")

pool = create_workers(verbose = true)
@everywhere begin
    using Arblib, CGL
    setprecision(Arb, 128)

    # Set logging to always flush
    using Logging: global_logger
    using TerminalLoggers: TerminalLogger
    global_logger(TerminalLogger(always_flush = true))
end

j, d, start, stop = read_args()
verbose = true

verbose && @info "Determined arguments" j d start stop

# TODO: Setup logging

verbose && @info "Computing initial point on branch"

μ, γ, κ, ξ₁, λ = CGL.sverak_params(Arb, j, d)

verbose && @info "Computing the branch"

br =
    let λ = CGL.CGLBranch.Params(λ.d, Float64(λ.ω), Float64(λ.σ), Float64(λ.δ), Float64(ξ₁))
        CGL.CGLBranch.branch_epsilon(Float64(μ), Float64(κ), Float64(λ.ϵ), λ)
    end

verbose && @info "Number of branch points" length(br)

if isnothing(stop)
    stop = length(br)
end

verbose && @info "Verifying branch"

success, ϵs, exists, uniqs = CGL.verify_branch(
    Arb.(br.param)[start:stop],
    Arb.(br.μ)[start:stop],
    Arb.(br.κ)[start:stop],
    ξ₁,
    λ,
    verbose = true,
    verbose_segments = true,
    log_progress = true,
)

verbose && @info "Finished verification of branch" length(ϵs) sum(success)
