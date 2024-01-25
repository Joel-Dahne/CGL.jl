"""
    classify_branch_parts(ϵs, μs, κs; cutoff = 10)

Return indices `i1` and `i2` such that the part where the derivative
of `κ` w.r.t. `ϵ` is greater than `cutoff` is contained between these
two indices.

With derivative we here mean the finite difference
```
(κs[i + 1] - κs[i]) / (ϵs[i + 1] - ϵs[i])
```
"""
function classify_branch_parts(ϵs::Vector{T}, κs::Vector{T}; cutoff::T = T(5)) where {T}
    dκ_dϵ = i -> (κs[i+1] - κs[i]) / (ϵs[i+1] - ϵs[i])

    i1 = findfirst(i -> abs(dκ_dϵ(i)) > cutoff, eachindex(ϵs, κs)[1:end-1])

    i2 = findlast(i -> abs(dκ_dϵ(i)) > cutoff, eachindex(ϵs, κs)[1:end-1])

    @assert !isnothing(i1) && !isnothing(i2)

    return i1, i2
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
            verbose = verbose_segments;
            maxevals,
            depth,
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
