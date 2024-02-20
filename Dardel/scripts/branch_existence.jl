using Dates

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
j, d, part, N = read_args()

verbose && @info "Determined arguments" j d part N

fix_kappa = part == "turn"

verbose && @info "Computing branch"

μ, γ, κ, ϵ, ξ₁, λ = CGL.sverak_params(Arb, j, d)

# We always want to use ξ₁ = 30 here
br = let λ = CGL.CGLBranch.Params(λ.d, λ.ω, λ.σ, λ.δ, 30.0)
    CGL.CGLBranch.branch_epsilon(Float64(μ), Float64(κ), Float64(ϵ), λ)
end

# Old version for choosing where to split between top, turn and bottom
#cutoff = ifelse(d == 3, 1.0, 10.0) # TODO: Tune this
#start_turning, stop_turning = CGL.classify_branch_parts(br.κ, br.param; cutoff)

start_turning, stop_turning = CGL.classify_branch_parts_2(br.param)

verbose && @info "Got $(length(br)) branch points" start_turning stop_turning

if part == "top"
    start = 1
    stop = start_turning
elseif part == "turn"
    # We want one overlapping segment with the top and bottom parts
    start = max(start_turning - 1, 1)
    stop = min(stop_turning + 1, length(br))

    if start == stop
        verbose && @error "Trying to compute turning part, but branch has no turning part"
    end
elseif part == "bottom"
    start = stop_turning
    stop = length(br)

    if start == stop
        verbose && @error "Trying to compute bottom part, but branch has no bottom part"
    end
else
    throw(ArgumentError("unknown part type $part"))
end

verbose && @info "Part \"$part\" has $(stop - start) segments"

stop_max = something(start_turning, length(br))

if !isnothing(N) && N < stop - start
    verbose && @info "Limiting to $N segments"
    stop = start + N
end

verbose && @info "Verifying branch"

runtime = @elapsed ϵs_or_κs, exists, uniqs, approxs = CGL.branch_existence(
    Arb.(br.μ[start:stop]),
    Arb.(br.κ[start:stop]),
    Arb.(br.param[start:stop]),
    ξ₁,
    λ,
    maxevals = 50000,
    verbose = true,
    verbose_segments = true;
    fix_kappa,
    pool,
)

dirname = "Dardel/output/branch_existence/$(round(Dates.now(), Second))"
filename = "branch_existence_j=$(j)_d=$(d)_part=$part.csv.gz"

verbose && @info "Writing data" dirname filename

if !fix_kappa
    df = CGL.branch_existence_dataframe_fix_epsilon(ϵs_or_κs, uniqs, exists, approxs)
else
    df = CGL.branch_existence_dataframe_fix_kappa(ϵs_or_κs, uniqs, exists, approxs)
end

mkpath(dirname)
CGL.write_parameters(joinpath(dirname, "parameters.csv"), ξ₁, λ; runtime)
CGL.write_branch_existence_csv(joinpath(dirname, filename), df)
