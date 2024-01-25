function verify_branch_existence(
    ϵs::Vector{Arf},
    μs::Vector{Arb},
    κs::Vector{Arb},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    pool = Distributed.WorkerPool(Distributed.workers()),
    verbose = false,
    verbose_segments = false,
    log_progress = false,
)
    @assert length(ϵs) == length(μs) == length(κs)

    verbose && @info "Verifying $(length(ϵs)-1) segments"

    tasks = map(1:length(ϵs)-1) do i
        @async Distributed.remotecall_fetch(
            verify_branch_segment_existence,
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

    failed_segments = map(segments) do (_, exists, _)
        !all(x -> all(isfinite, x), exists)
    end

    if iszero(any(failed_segments))
        verbose && @info "Succesfully verified all segments"
    else
        verbose && @warn "Failed verifying $(sum(failed_segments)) segments"
    end

    ϵs = reduce(vcat, getindex.(segments, 1))
    exists = reduce(vcat, getindex.(segments, 2))
    uniqs = reduce(vcat, getindex.(segments, 3))

    @assert issorted(ϵs, by = x -> x[1])

    return ϵs, exists, uniqs
end
