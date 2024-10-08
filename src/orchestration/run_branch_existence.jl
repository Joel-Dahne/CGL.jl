function run_branch_existence(
    j::Integer = 1,
    d::Integer = 1,
    part = "top";
    try_expand_uniqueness::Bool = true,
    indices::Union{Nothing,UnitRange{Int}} = nothing,
    N::Union{Nothing,Integer} = nothing,
    maxevals::Integer = 50000,
    depth::Integer = 20,
    pool::Distributed.WorkerPool = Distributed.WorkerPool(Distributed.workers()),
    save_results::Bool = true,
    directory::Union{Nothing,AbstractString} = nothing,
    log_progress::Bool = true,
    verbose::Bool = true,
    verbose_segments::Bool = true,
)
    verbose && @info "Computing for j = $j,  d = $d, part = $part" N

    fix_kappa = part == "turn"

    verbose && @info "Computing initial branch"

    _, _, _, _, ξ₁, λ = CGL.sverak_params(Arb, j, d, ξ₁_for_branch = true)

    μs, κs, ϵs = let
        br = CGL.CGLBranch.branch_epsilon(CGL.CGLBranch.sverak_initial(j, d)...)

        Arb.(br.μ), Arb.(br.κ), Arb.(br.param)
    end

    if isnothing(indices)
        start_turning, stop_turning = CGL.classify_branch_parts(ϵs)

        verbose && @info "Got $(length(μs)) branch points" start_turning stop_turning

        if part == "top"
            start = 1
            stop = start_turning
        elseif part == "turn"
            # We want one overlapping segment with the top and bottom parts
            start = max(start_turning - 1, 1)
            stop = min(stop_turning + 1, length(μs))

            if start == stop
                verbose &&
                    @error "Trying to compute turning part, but branch has no turning part"
            end
        elseif part == "bottom"
            start = stop_turning
            stop = length(μs)

            if start == stop
                verbose &&
                    @error "Trying to compute bottom part, but branch has no bottom part"
            end
        else
            throw(ArgumentError("unknown part $part"))
        end

        verbose && @info "Part \"$part\" has $(stop - start) segments"

        stop_max = something(start_turning, length(μs))

        if !isnothing(N) && N < stop - start
            verbose && @info "Limiting to $N segments"
            stop = start + N
        end

        indices = start:stop
    end

    verbose && @info "Verifying branch"

    runtime = @elapsed ϵs_or_κs, exists, uniqs, approxs = CGL.branch_existence(
        μs[indices],
        κs[indices],
        ϵs[indices],
        ξ₁,
        λ;
        fix_kappa,
        pool,
        maxevals,
        depth,
        try_expand_uniqueness,
        verbose,
        verbose_segments,
        log_progress,
    )

    if !fix_kappa
        df = CGL.branch_existence_dataframe_fix_epsilon(ϵs_or_κs, uniqs, exists, approxs)
    else
        df = CGL.branch_existence_dataframe_fix_kappa(ϵs_or_κs, uniqs, exists, approxs)
    end

    if save_results
        if isnothing(directory)
            # Avoid colon (:) in date string since that is not allowed
            # on Windows
            date_string = replace(string(round(Dates.now(), Dates.Second)), ":" => "")
            commit_string = readchomp(`git rev-parse --short HEAD`)
            directory = relpath(
                joinpath(
                    dirname(pathof(@__MODULE__)),
                    "../HPC/output/branch_existence_d=$(d)_j=$(j)_part=$(part)",
                    "$(date_string)_$(commit_string)",
                ),
            )
        end

        verbose && @info "Writing data" directory

        mkpath(directory)
        CGL.write_parameters(
            joinpath(directory, "parameters.csv"),
            ξ₁,
            λ;
            runtime,
            commit_hash = readchomp(`git rev-parse HEAD`),
        )
        CGL.write_branch_csv(
            joinpath(directory, "branch_existence_d=$(d)_j=$(j)_part=$part.csv.gz"),
            df,
        )
    else
        verbose && @info "Not writing data"
    end

    return df
end
