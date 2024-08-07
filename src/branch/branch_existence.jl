"""
    classify_branch_parts(ϵs)

Return indices splitting the vector into a top part, a turning part
and a bottom part.

It return two indices, `start_turning` and `stop_turning`. The top
part of the branch is given by the indices `1:start_turning`, the turn
part by `start_turning:stop_turning` and the bottom part by
`stop_turning:lastindex(ϵs)`.

And index range contain only one index means that this part was not
present on the branch (for example because it stops before the turning
point).

The indices are determined by first finding the `ϵ` value right before
the turn. The top part is then all indices for which the `ϵ` value is
less than half that at the turn and the bottom part is all indices
where the `ϵ` value is halfway between the turn and the endpoint.
"""
function classify_branch_parts(ϵs::Vector{T}) where {T}
    turning_point = findfirst(i -> ϵs[i+1] < ϵs[i], eachindex(ϵs)[1:end-1])

    # Branch doesn't reach the turning point. Classify all of it as
    # belonging to the top part.
    isnothing(turning_point) && return length(ϵs), length(ϵs)

    ϵ_turning_point = ϵs[turning_point]

    start_turning = findfirst(i -> ϵs[i] > ϵ_turning_point / 2, eachindex(ϵs))
    isnothing(start_turning) && error("could not determinate start of turning")

    stop_turning = findlast(i -> ϵs[i] > (ϵ_turning_point + ϵs[end]) / 2, eachindex(ϵs))
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
    try_expand_uniqueness = true,
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
                try_expand_uniqueness,
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
                try_expand_uniqueness,
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
