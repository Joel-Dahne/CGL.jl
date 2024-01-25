using CSV, DataFrames, Dates

include("helper.jl")

pool, num_threads = create_workers(verbose = true)

@everywhere begin
    using Arblib, CGL
    setprecision(Arb, 128)

    # Set logging to always flush
    using Logging: global_logger
    using TerminalLoggers: TerminalLogger
    global_logger(TerminalLogger(always_flush = true))
end

verbose = true

# Read arguments
j, d, start, stop = read_args()

verbose && @info "Determined arguments" j d start stop

verbose && @info "Computing branch"

μ, γ, κ, ξ₁, λ = CGL.sverak_params(Arb, j, d)

# We always want to use ξ₁ = 30 here
br = let λ = CGL.CGLBranch.Params(λ.ϵ, λ.d, λ.ω, λ.σ, λ.δ, 30.0)
    CGL.CGLBranch.branch(Float64(μ), Float64(κ), λ)
end

verbose && @info "Got $(length(br.branch)) branch points"

start_turning, _ = CGL.classify_branch_parts(br.branch.param, br.branch.κ)

stop_max = something(start_turning, length(br.branch))

verbose && @info "$stop_max segments before turning starts"

if !isnothing(stop) && stop < stop_max
    verbose && @info "Limiting to $(stop - start) segments"
else
    stop = stop_max
end

verbose && @info "Verifying branch"

ϵs, exists, uniqs = CGL.verify_branch_existence(
    Arf.(br.branch.param)[start:stop],
    Arb.(br.branch.μ)[start:stop],
    Arb.(br.branch.κ)[start:stop],
    ξ₁,
    λ,
    verbose = true,
    verbose_segments = true,
    log_progress = true;
    pool,
)

dirname = "Dardel/output/branch_existence/$(round(Dates.now(), Second))"
mkpath(dirname)

verbose && @info "Writing data" dirname

df = CGL.branch_dataframe(
    ϵs,
    getindex.(uniqs, 1),
    Acb.(getindex.(uniqs, 2), getindex.(uniqs, 3)),
    getindex.(uniqs, 4),
    getindex.(exists, 1),
    Acb.(getindex.(exists, 2), getindex.(exists, 3)),
    getindex.(exists, 4),
    fill(ξ₁, length(ϵs)),
)

filename = "branch_existence_j=$(j)_d=$d.csv"
CGL.write_branch_csv(joinpath(dirname, filename), df)
