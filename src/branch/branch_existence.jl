"""
    classify_branch_parts(ϵs, μs, κs; cutoff = 10)

Return indices `start_turning` and `stop_turning` such that the
turning point is contained between these two indices and such that the
derivative of `κ` w.r.t. `ϵ` is greater than `cutoff` between them.

With derivative we here mean the finite difference
```
(κs[i + 1] - κs[i]) / (ϵs[i + 1] - ϵs[i])
```
"""
function classify_branch_parts(ϵs::Vector{T}, κs::Vector{T}; cutoff::T = T(1)) where {T}
    dκ_dϵ = i -> (κs[i+1] - κs[i]) / (ϵs[i+1] - ϵs[i])

    turning_point = findfirst(i -> dκ_dϵ(i) > 0, eachindex(ϵs, κs)[1:end-1])

    if isnothing(turning_point)
        # Branch doesn't reach the turning point
        start_turning = findlast(i -> abs(dκ_dϵ(i)) < cutoff, eachindex(ϵs, κs)[1:end-1])

        if isnothing(start_turning)
            # Branch doesn't get close to turning point
            return length(ϵs), length(ϵs)
        else
            return start_turning, length(ϵs)
        end
    end

    if turning_point == 1
        # If this happens it most likely means that we were not able
        # to accurate compute the branch.
        error("branch is not decreasing at the start")
    end

    start_turning = findlast(i -> abs(dκ_dϵ(i)) < cutoff, 1:turning_point-1)

    # There should be at least some segments around the turning point
    @assert !isnothing(start_turning) && 1 < start_turning < turning_point - 1

    stop_turning_part = findfirst(i -> abs(dκ_dϵ(i)) < cutoff, turning_point:length(ϵs)-1)

    if isnothing(stop_turning_part)
        stop_turning = length(ϵs)
    else
        stop_turning = turning_point - 1 + stop_turning_part
    end

    return start_turning, stop_turning
end

function verify_branch_existence(
    ϵs::Vector{Arf},
    μs::Vector{Arb},
    κs::Vector{Arb},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    pool = Distributed.WorkerPool(Distributed.workers()),
    maxevals::Integer = 1000,
    depth::Integer = 20,
    verbose = false,
    verbose_segments = false,
    log_progress = verbose,
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
            verbose = verbose_segments;
            maxevals,
            depth,
        )
    end

    segments = fetch_with_progress(tasks, log_progress)

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
    approxs = reduce(vcat, getindex.(segments, 4))

    return ϵs, exists, uniqs, approxs
end

function verify_branch_existence_epsilon(
    ϵs::Vector{Arb},
    μs::Vector{Arb},
    κs::Vector{Arf},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    pool = Distributed.WorkerPool(Distributed.workers()),
    maxevals::Integer = 1000,
    depth::Integer = 20,
    verbose = false,
    verbose_segments = false,
    log_progress = verbose,
)
    @assert length(ϵs) == length(μs) == length(κs)

    verbose && @info "Verifying $(length(κs)-1) segments"

    tasks = map(1:length(κs)-1) do i
        @async Distributed.remotecall_fetch(
            verify_branch_segment_existence_epsilon,
            pool,
            (ϵs[i], ϵs[i+1]),
            (μs[i], μs[i+1]),
            (κs[i], κs[i+1]),
            ξ₁,
            λ,
            verbose = verbose_segments;
            maxevals,
            depth,
        )
    end

    segments = fetch_with_progress(tasks, log_progress)

    failed_segments = map(segments) do (_, exists, _)
        !all(x -> all(isfinite, x), exists)
    end

    if iszero(any(failed_segments))
        verbose && @info "Succesfully verified all segments"
    else
        verbose && @warn "Failed verifying $(sum(failed_segments)) segments"
    end

    κs = reduce(vcat, getindex.(segments, 1))
    exists = reduce(vcat, getindex.(segments, 2))
    uniqs = reduce(vcat, getindex.(segments, 3))
    approxs = reduce(vcat, getindex.(segments, 4))

    return κs, exists, uniqs, approxs
end
