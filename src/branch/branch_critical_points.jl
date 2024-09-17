function _branch_critical_points_batch_mince(
    μ::Arb,
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    fix_kappa::Bool = false,
    max_depth::Integer = 0,
)
    for depth = 1:max_depth
        minced = if fix_kappa
            mince(κ, 2^depth)
        else
            mince(ϵ, 2^depth)
        end

        num_critical_points =
            convert(Vector{Union{Missing,Integer}}, fill(missing, length(minced)))
        @show depth
        for j in eachindex(minced)
            if fix_kappa
                κ_minced = minced[j]
                G_x = x -> G(x[1:3]..., κ_minced, x[4], ξ₁, λ)
                dG_x = x -> G_jacobian_epsilon(x[1:3]..., κ_minced, x[4], ξ₁, λ)

                μ_minced, γ_real_minced, γ_imag_minced, ϵ_minced =
                    verify_and_refine_root(G_x, dG_x, SVector(µ, real(γ), imag(γ), ϵ))
            else
                ϵ_minced = minced[j]
                G_x = x -> G(x..., ϵ_minced, ξ₁, λ)
                dG_x = x -> G_jacobian_kappa(x..., ϵ_minced, ξ₁, λ)

                μ_minced, γ_real_minced, γ_imag_minced, κ_minced =
                    verify_and_refine_root(G_x, dG_x, SVector(µ, real(γ), imag(γ), κ))
            end

            success, zeros, verified_zeros = count_critical_points(
                μ_minced,
                Acb(γ_real_minced, γ_imag_minced),
                κ_minced,
                ϵ_minced,
                ξ₁,
                λ,
            )

            if success
                num_critical_points[j] = length(verified_zeros)
            else
                break
            end
        end

        if all(!ismissing, num_critical_points)
            return only(unique(num_critical_points))
        end
    end

    # If we got here it didn't succeed for any used depth
    return missing
end

function branch_critical_points_batch(
    μs::Vector{Arb},
    γs::Vector{Acb},
    κs::Vector{Arb},
    ϵs::Vector{Arb},
    ξ₁s::Vector{Arb},
    λ::CGLParams{Arb};
    fix_kappas::Vector{Bool} = fill(false, length(μs)),
    max_depth::Integer = 6,
    verbose = false,
)
    res = tmap(
        Union{Int,Missing},
        eachindex(μs, γs, κs, ϵs, ξ₁s, fix_kappas),
        scheduler = :greedy,
    ) do i
        success, zeros, verified_zeros =
            count_critical_points(μs[i], γs[i], κs[i], ϵs[i], ξ₁s[i], λ)

        if success
            length(verified_zeros)
        elseif max_depth == 0
            missing
        else
            _branch_critical_points_batch_mince(
                μs[i],
                γs[i],
                κs[i],
                ϵs[i],
                ξ₁s[i],
                λ;
                fix_kappa = fix_kappas[i],
                max_depth,
            )
        end
    end

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
    fix_kappas::Vector{Bool} = fill(false, length(μs)),
    max_depth::Integer = 6,
    pool = Distributed.WorkerPool(Distributed.workers()),
    batch_size = 128,
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
            fix_kappas = fix_kappas[indices_batch],
            max_depth,
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
                @info "Got $(res_ok[1]) critical points for all $(ifelse(num_missing  > 0, "successfull ", ""))parts"
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
