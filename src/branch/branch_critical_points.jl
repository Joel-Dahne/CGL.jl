function branch_critical_points_batch(
    μs::Vector{Arb},
    γs::Vector{Acb},
    κs::Vector{Arb},
    ϵs::Vector{Arb},
    ξ₁s::Vector{Arb},
    λ::CGLParams{Arb};
    use_midpoint::Bool = false,
    verbose = false,
)
    res = tmap(Union{Int,Missing}, eachindex(μs, γs, κs, ϵs), scheduler = :greedy) do i
        success, zeros, verified_zeros =
            count_critical_points(μs[i], γs[i], κs[i], ϵs[i], ξ₁s[i], λ; use_midpoint)

        ifelse(success, length(verified_zeros), missing)
    end

    # There has been issues with high memory consumption giving
    # OOM crashes on SLURM. Explicitly galling gc here helps with
    # that.
    GC.gc()

    if verbose
        if any(ismissing, res)
            @warn "Couldn't compute number of critical points for full batch"
        elseif allequal(res)
            @info "Got $(res[1]) critical points for full batch"
        else
            # This shouldn't happen!
            @error "Not all parts of patch got same number of critical points" res
        end
    end

    return res
end

function branch_critical_points(
    μs::Vector{Arb},
    γs::Vector{Acb},
    κs::Vector{Arb},
    ϵs::Vector{Arb},
    ξ₁s::Vector{Arb},
    λ::CGLParams{Arb};
    use_midpoint::Bool = false,
    pool = Distributed.WorkerPool(Distributed.workers()),
    batch_size = 32,
    verbose = false,
    log_progress = verbose,
)
    # If the batch_size is too large then not all workers will get a
    # share. Take it smaller if this is the case.
    batch_size = min(batch_size, ceil(Int, length(μs) / length(pool)))

    indices = firstindex(μs):batch_size:lastindex(μs)

    verbose && @info "Starting $(length(indices)) batch jobs of size $batch_size"

    tasks = map(indices) do index
        indices_batch = index:min(index + batch_size - 1, lastindex(μs))

        @async Distributed.remotecall_fetch(
            (args...; kwargs...) ->
                @time(branch_critical_points_batch(args...; kwargs...)),
            pool,
            μs[indices_batch],
            γs[indices_batch],
            κs[indices_batch],
            ϵs[indices_batch],
            ξ₁s[indices_batch],
            λ;
            use_midpoint,
            verbose,
        )
    end

    verbose && @info "Collecting batch jobs"

    batches = fetch_with_progress(tasks, log_progress)

    res = foldl(vcat, batches)

    if verbose
        num_missing = count(ismissing, res)

        if num_missing > 0
            @warn "Failed computing number of critical points for $num_missing / $(length(res)) parts"
        end

        res_ok = filter(!ismissing, res)

        if !isempty(res_ok)
            if allequal(res_ok)
                @info "Got $(res[1]) critical points for all $(ifelse(num_missing  > 0, "successfull ", ""))parts"
            else
                # This shouldn't happen!
                @error "Not all parts  got same number of critical points" unique(res) findall(
                    !=(res[1]),
                    res,
                )
            end
        end
    end

    return res
end
