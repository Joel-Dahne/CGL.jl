function verify_branch_points_batch(
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

function verify_branch_points_batch_epsilon(
    μs::AbstractVector{Arb},
    κs::AbstractVector{Arb},
    ξ₁s::AbstractVector{Arb},
    λs::AbstractVector{CGLParams{Arb}},
)
    res = similar(μs, SVector{4,Arb})

    Threads.@threads :dynamic for i in eachindex(res, μs, κs, ξ₁s, λs)
        μ, γ, ϵ = refine_approximation_epsilon(μs[i], κs[i], ξ₁s[i], λs[i])

        if abs(μ - μs[i]) > 0.1 || abs(ϵ - λs[i].ϵ) > 0.1
            # If we are this far off from the initial approximation
            # then something went wrong.
            res[i] = SVector(
                indeterminate(μ),
                indeterminate(μ),
                indeterminate(μ),
                indeterminate(μ),
            )
        else
            res[i] =
                G_solve_epsilon(μ, real(γ), imag(γ), κs[i], ξ₁s[i], CGLParams(λs[i]; ϵ))
        end
    end

    return res
end

function verify_branch_points(
    μs::Vector{Arb},
    κs::Vector{Arb},
    ξ₁s::Vector{Arb},
    λs::Vector{CGLParams{Arb}};
    use_epsilon = false,
    pool = Distributed.WorkerPool(Distributed.workers()),
    batch_size = 128,
    verbose = false,
    log_progress = false,
)
    @assert length(μs) == length(κs) == length(ξ₁s) == length(λs)

    indices = firstindex(μs):batch_size:lastindex(μs)

    verbose && @info "Starting $(length(indices)) batch jobs of size $batch_size"

    tasks = map(indices) do index
        indices_batch = index:min(index + batch_size - 1, lastindex(μs))

        if !use_epsilon
            @async Distributed.remotecall_fetch(
                verify_branch_points_batch,
                pool,
                μs[indices_batch],
                κs[indices_batch],
                ξ₁s[indices_batch],
                λs[indices_batch],
            )
        else
            @async Distributed.remotecall_fetch(
                verify_branch_points_batch_epsilon,
                pool,
                μs[indices_batch],
                κs[indices_batch],
                ξ₁s[indices_batch],
                λs[indices_batch],
            )
        end
    end

    verbose && @info "Collecting batch jobs"

    batches = fetch_with_progress(tasks, log_progress)

    res = foldl(vcat, batches)

    return res
end
