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

verbose && @info "Computing branch"

μ, γ, κ, ξ₁, λ = CGL.sverak_params(Arb, j, d)

# We always want to use ξ₁ = 30 here
br = let ϵ = λ.ϵ, λ = CGL.CGLBranch.Params(λ.d, λ.ω, λ.σ, λ.δ, 30.0)
    CGL.CGLBranch.branch_epsilon(Float64(μ), Float64(κ), Float64(ϵ), λ)
end

cutoff = ifelse(d == 3, 1.0, 10.0)

start_turning, stop_turning = CGL.classify_branch_parts(br.param, br.κ; cutoff)

verbose && @info "Got $(length(br)) branch points" start_turning stop_turning

if part == "top"
    start = 1
    stop = start_turning
elseif part == "turn"
    start = start_turning
    stop = stop_turning

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

if part == "top" || part == "bottom"
    ϵs, exists, uniqs, approxs = CGL.verify_branch_existence(
        Arf.(br.param)[start:stop],
        Arb.(br.μ)[start:stop],
        Arb.(br.κ)[start:stop],
        ξ₁,
        λ,
        maxevals = 5000,
        verbose = true,
        verbose_segments = true;
        pool,
    )
elseif part == "turn"
    κs, exists, uniqs, approxs = CGL.verify_branch_existence_epsilon(
        Arb.(br.param)[start:stop],
        Arb.(br.μ)[start:stop],
        Arf.(br.κ)[start:stop],
        ξ₁,
        λ,
        maxevals = 5000,
        verbose = true,
        verbose_segments = true;
        pool,
    )
else
    throw(ArgumentError("unknown part type $part"))
end

dirname = "Dardel/output/branch_existence/$(round(Dates.now(), Second))"
filename = "branch_existence_j=$(j)_d=$(d)_part=$part.csv"
mkpath(dirname)

verbose && @info "Writing data" dirname filename

if part == "top" || part == "bottom"
    df = CGL.branch_existence_dataframe(ϵs, uniqs, exists, approxs, ξ₁)
elseif part == "turn"
    df = CGL.branch_existence_dataframe_epsilon(κs, uniqs, exists, approxs, ξ₁)
else
    throw(ArgumentError("unknown part type $part"))
end

CGL.write_branch_existence_csv(joinpath(dirname, filename), df)
