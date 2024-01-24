function verify_branch(
    ϵs::Vector{Arb},
    μs::Vector{Arb},
    κs::Vector{Arb},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    verbose = false,
    verbose_segments = false,
    log_progress = false,
)
    @assert length(ϵs) == length(μs) == length(κs)

    pool = Distributed.WorkerPool(Distributed.workers())

    verbose && @info "Verifying $(length(ϵs)-1) segments"

    tasks = map(1:length(ϵs)-1) do i
        @async Distributed.remotecall_fetch(
            verify_branch_segment,
            pool,
            (ϵs[i], ϵs[i+1]),
            (μs[i], μs[i+1]),
            (κs[i], κs[i+1]),
            ξ₁,
            λ,
            verbose = verbose_segments,
        )
    end

    if log_progress
        @progress segments = [fetch(task) for task in tasks]
    else
        segments = [fetch(task) for task in tasks]
    end

    failed_segments = map(segments) do (_, exists, uniqs)
        any(check_continuation(exists, uniqs))
    end

    if iszero(any(failed_segments))
        verbose && @info "Succesfully verified all segments"
    else
        verbose && @warn "Failed verifying $(sum(failed_segments)) segments"
    end

    ϵs = reduce(vcat, getindex.(segments, 1))
    exists = reduce(vcat, getindex.(segments, 2))
    uniqs = reduce(vcat, getindex.(segments, 3))

    to_bisect = check_continuation(exists, uniqs)

    if iszero(any(failed_segments))
        # Check if any endpoints needs to be bisected
        if any(to_bisect)
            verbose &&
                @info "Continuation failed for $(sum(to_bisect)) intersections of segments"

            # TODO: Implement this
            verbose && @error "Bisection for intersections not yet implemented"
        else
            verbose && @info "Verified continuation for all intersections of segments"
        end
    else
        verbose && @warn "Not attempting to veryify intersection of segments"
    end

    success = .!check_continuation(exists, uniqs)

    return success, ϵs, exists, uniqs
end
