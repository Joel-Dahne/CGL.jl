function branch_points_batch_fix_epsilon(
    μs::AbstractVector{Arb},
    κs::AbstractVector{Arb},
    ϵs::AbstractVector{Arb},
    ξ₁s::AbstractVector{<:Union{Arb,AbstractVector{Arb}}},
    λs::AbstractVector{CGLParams{Arb}},
)
    return tmap(
        Tuple{SVector{4,Arb},Arb},
        eachindex(μs, κs, ϵs, ξ₁s, λs),
        scheduler = :greedy,
    ) do i
        for ξ₁ in ξ₁s[i] # Abuse that a number can iterated as a singleton vector
            converged, μ, γ, κ = refine_approximation_fix_epsilon(
                μs[i],
                κs[i],
                ϵs[i],
                ξ₁,
                λs[i],
                return_convergence = Val{true}(),
            )

            # Avoid calling G_solve the refinement didn't converge or
            # if we are very far from the initial approximation.
            # Something probably went wrong in that case.
            if converged && abs(μ - μs[i]) < 0.1 && abs(κ - κs[i]) < 0.1
                exist = G_solve_fix_epsilon_alt(μ, real(γ), imag(γ), κ, ϵs[i], ξ₁, λs[i])

                all(isfinite, exist) && return exist, ξ₁
            end
        end

        # Return an indeterminate value if no ξ₁ worked
        exist = SVector(
            indeterminate(Arb),
            indeterminate(Arb),
            indeterminate(Arb),
            indeterminate(Arb),
        )
        return exist, indeterminate(Arb)
    end
end

function branch_points_batch_fix_kappa(
    μs::AbstractVector{Arb},
    κs::AbstractVector{Arb},
    ϵs::AbstractVector{Arb},
    ξ₁s::AbstractVector{<:Union{Arb,AbstractVector{Arb}}},
    λs::AbstractVector{CGLParams{Arb}},
)
    return tmap(
        Tuple{SVector{4,Arb},Arb},
        eachindex(μs, κs, ϵs, ξ₁s, λs),
        scheduler = :greedy,
    ) do i
        for ξ₁ in ξ₁s[i] # Abuse that a number can iterated as a singleton vector
            converged, μ, γ, ϵ = refine_approximation_fix_kappa(
                μs[i],
                κs[i],
                ϵs[i],
                ξ₁,
                λs[i],
                return_convergence = Val{true}(),
            )

            # Avoid calling G_solve the refinement didn't converge or
            # if we are very far from the initial approximation.
            # Something probably went wrong in that case.
            if converged && abs(μ - μs[i]) < 0.1 && abs(ϵ - ϵs[i]) < 0.1
                exist = G_solve_fix_kappa_alt(μ, real(γ), imag(γ), κs[i], ϵ, ξ₁, λs[i])

                all(isfinite, exist) && return exist, ξ₁
            end
        end

        # Return an indeterminate value if no ξ₁ worked
        exist = SVector(
            indeterminate(Arb),
            indeterminate(Arb),
            indeterminate(Arb),
            indeterminate(Arb),
        )
        return exist, indeterminate(Arb)
    end
end

function branch_points(
    μs::Vector{Arb},
    κs::Vector{Arb},
    ϵs::Vector{Arb},
    ξ₁s::Vector{<:Union{Arb,AbstractVector{Arb}}},
    λs::Vector{CGLParams{Arb}};
    fix_kappa = false,
    pool = Distributed.WorkerPool(Distributed.workers()),
    batch_size = 128,
    verbose = false,
    log_progress = false,
)
    @assert length(μs) == length(κs) == length(ϵs) == length(ξ₁s) == length(λs)

    indices = firstindex(μs):batch_size:lastindex(μs)

    verbose && @info "Starting $(length(indices)) batch jobs of size $batch_size"

    tasks = map(indices) do index
        indices_batch = index:min(index + batch_size - 1, lastindex(μs))

        if !fix_kappa
            @async Distributed.remotecall_fetch(
                (args...; kwargs...) ->
                    @time(branch_points_batch_fix_epsilon(args...; kwargs...)),
                pool,
                μs[indices_batch],
                κs[indices_batch],
                ϵs[indices_batch],
                ξ₁s[indices_batch],
                λs[indices_batch],
            )
        else
            @async Distributed.remotecall_fetch(
                (args...; kwargs...) ->
                    @time(branch_points_batch_fix_kappa(args...; kwargs...)),
                pool,
                μs[indices_batch],
                κs[indices_batch],
                ϵs[indices_batch],
                ξ₁s[indices_batch],
                λs[indices_batch],
            )
        end
    end

    verbose && @info "Collecting batch jobs"
    batches = fetch_with_progress(tasks, log_progress)

    res = foldl(vcat, batches)

    exists = getindex.(res, 1)
    ξ₁s_used = getindex.(res, 2)

    if verbose
        success = count(isfinite.(getindex.(exists, 1)))
        @info "Succesfully verified $success / $(length(exists)) points"
    end

    return exists, ξ₁s_used
end
