"""
    classify_branch_parts(κs, ϵs; cutoff = 10)

Return indices `start_turning` and `stop_turning` such that the
turning point is contained between these two indices and such that the
derivative of `κ` w.r.t. `ϵ` is greater than `cutoff` between them.

With derivative we here mean the finite difference
```
(κs[i + 1] - κs[i]) / (ϵs[i + 1] - ϵs[i])
```
"""
function classify_branch_parts(κs::Vector{T}, ϵs::Vector{T}; cutoff::T = T(1)) where {T}
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

"""
    classify_branch_parts_2(ϵs)

Alternative to [`classify_branch_parts`](@ref) which instead returns
the midpoint of the top and bottom part of the branch respectively.
"""
function classify_branch_parts_2(ϵs::Vector{T}) where {T}
    dκ_dϵ = i -> (κs[i+1] - κs[i]) / (ϵs[i+1] - ϵs[i])

    turning_point = findfirst(i -> ϵs[i+1] < ϵs[i], eachindex(κs)[1:end-1])

    if isnothing(turning_point)
        # Branch doesn't reach the turning point. Classify all of it
        # as belonging to the top part.

        return length(ϵs), length(ϵs)
    end

    start_turning = findfirst(i -> ϵs[i] > ϵs[turning_point] / 2, eachindex(ϵs))
    isnothing(start_turning) && error("could not determinate start of turning")

    stop_turning = findlast(i -> ϵs[i] > (ϵs[turning_point] + ϵs[end]) / 2, eachindex(ϵs))
    isnothing(stop_turning) && error("could not determinate stop of turning")

    return start_turning, stop_turning
end

function branch_existence(
    μs::Vector{Arb},
    κs::Vector{Arb},
    ϵs::Vector{Arb},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    fix_kappa = false,
    pool = Distributed.WorkerPool(Distributed.workers()),
    maxevals::Integer = 1000,
    depth::Integer = 20,
    verbose = false,
    verbose_segments = false,
    log_progress = verbose,
)
    @assert length(μs) == length(κs) == length(ϵs)

    verbose && @info "Verifying $(length(μs)-1) segments"

    tasks = map(1:length(μs)-1) do i
        if !fix_kappa
            @async Distributed.remotecall_fetch(
                (args...; kwargs...) ->
                    @time(branch_segment_existence_fix_epsilon(args...; kwargs...)),
                pool,
                (μs[i], μs[i+1]),
                (κs[i], κs[i+1]),
                (midpoint(ϵs[i]), midpoint(ϵs[i+1])),
                ξ₁,
                λ,
                verbose = verbose_segments;
                maxevals,
                depth,
            )
        else
            @async Distributed.remotecall_fetch(
                (args...; kwargs...) ->
                    @time(branch_segment_existence_fix_kappa(args...; kwargs...)),
                pool,
                (μs[i], μs[i+1]),
                (midpoint(κs[i]), midpoint(κs[i+1])),
                (ϵs[i], ϵs[i+1]),
                ξ₁,
                λ,
                verbose = verbose_segments;
                maxevals,
                depth,
            )
        end
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

    ϵs_or_κs = reduce(vcat, getindex.(segments, 1))
    exists = reduce(vcat, getindex.(segments, 2))
    uniqs = reduce(vcat, getindex.(segments, 3))
    approxs = reduce(vcat, getindex.(segments, 4))

    return ϵs_or_κs, exists, uniqs, approxs
end
