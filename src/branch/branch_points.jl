function verify_branch_points(
    ϵs::Vector{Arb},
    μs::Vector{Arb},
    κs::Vector{Arb},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    verbose = false,
    log_progress = false,
)
    res = similar(ϵs, SVector{4,Arb})

    progress = Threads.Atomic{Int}(0)

    @withprogress Threads.@threads for i in eachindex(ϵs, μs, κs)
        λ_ϵ = CGLParams(λ, ϵ = ϵs[i])

        μ, γ, κ = refine_approximation(μs[i], κs[i], ξ₁, λ_ϵ)

        res[i] = G_solve(μ, real(γ), imag(γ), κ, ξ₁, λ_ϵ)

        Threads.atomic_add!(progress, 1)
        log_progress && @logprogress progress[] / length(res)
    end

    return res
end

function verify_branch_points_distributed_helper(
    μs::AbstractVector{Arb},
    κs::AbstractVector{Arb},
    ξ₁s::AbstractVector{Arb},
    λs::AbstractVector{CGLParams{Arb}},
)
    res = similar(μs, SVector{4,Arb})

    Threads.@threads :dynamic for i in eachindex(res, μs, κs, ξ₁s, λs)
        μ, γ, κ = refine_approximation(μs[i], κs[i], ξ₁s[i], λs[i])

        if abs(μ - μs[i]) > 0.1 || abs(κ - κs[i]) > 0.1
            # If we are this far off from the initial approximation
            # then something went wrong.
            res[i] = SVector(
                indeterminate(μ),
                indeterminate(μ),
                indeterminate(μ),
                indeterminate(μ),
            )
        else
            res[i] = G_solve(μ, real(γ), imag(γ), κ, ξ₁s[i], λs[i])
        end
    end

    return res
end

function verify_branch_points_distributed(
    μs::Vector{Arb},
    κs::Vector{Arb},
    ξ₁s::Vector{Arb},
    λs::Vector{CGLParams{Arb}};
    batch_size = 128,
    verbose = false,
    log_progress = false,
)
    @assert length(μs) == length(κs) == length(ξ₁s) == length(λs)

    pool = Distributed.WorkerPool(Distributed.workers())

    indices = firstindex(μs):batch_size:lastindex(μs)

    verbose && @info "Starting $(length(indices)) batch jobs of size $batch_size"

    tasks = map(indices) do index
        indices_batch = index:min(index + batch_size - 1, lastindex(μs))

        @async Distributed.remotecall_fetch(
            verify_branch_points_distributed_helper,
            pool,
            μs[indices_batch],
            κs[indices_batch],
            ξ₁s[indices_batch],
            λs[indices_batch],
        )
    end

    verbose && @info "Collecting batch jobs"

    if log_progress
        @progress batches = [fetch(task) for task in tasks]
    else
        batches = [fetch(task) for task in tasks]
    end

    res = foldl(vcat, batches)

    return res
end
