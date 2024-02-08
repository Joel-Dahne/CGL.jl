function check_left_continuation(
    exists::Vector{SVector{4,Arb}},
    uniqs::Vector{SVector{4,Arb}},
)
    left_continuation_finite = fill(false, length(exists))
    left_continuation = fill(false, length(exists))

    left_continuation_finite[1] = all(isfinite, exists[1])
    left_continuation[1] = all(isfinite, exists[1])

    for i in eachindex(exists, uniqs)[2:end]
        left_continuation_finite[i] = all(isfinite, uniqs[i-1]) && all(isfinite, exists[i])

        left_continuation[i] =
            left_continuation_finite[i] && all(Arblib.contains.(uniqs[i-1], exists[i]))
    end

    return left_continuation_finite, left_continuation
end

function continuation_insert_bisected(
    ϵs_or_κs::Vector{NTuple{2,Arf}},
    exists::Vector{SVector{4,Arb}},
    uniqs::Vector{SVector{4,Arb}},
    approxs::Vector{SVector{4,Arb}},
    to_bisect::Union{Vector{Bool},BitVector},
    ϵs_or_κs_bisected::Vector{NTuple{2,Arf}},
    exists_bisected::Vector{SVector{4,Arb}},
    uniqs_bisected::Vector{SVector{4,Arb}},
    approxs_bisected::Vector{SVector{4,Arb}},
)
    N = length(ϵs) + sum(to_bisect)

    ϵs_or_κs_new = typeof(ϵs)(undef, N)
    exists_new = typeof(exists)(undef, N)
    uniqs_new = typeof(uniqs)(undef, N)
    approxs_new = typeof(approxs)(undef, N)

    j = 1
    k = 1
    for i in eachindex(ϵs)
        if to_bisect[i]
            ϵs_or_κs_new[j:j+1] .= ϵs_or_κs_bisected[k:k+1]
            exists_new[j:j+1] .= exists_bisected[k:k+1]
            uniqs_new[j:j+1] .= uniqs_bisected[k:k+1]
            approxs_new[j:j+1] .= approxs_bisected[k:k+1]

            j += 2
            k += 2
        else
            ϵs_or_κs_new[j] = ϵs_or_κs[i]
            exists_new[j] = exists[i]
            uniqs_new[j] = uniqs[i]
            approxs_new[j] = approxs[i]

            j += 1
        end
    end

    @assert j == N + 1

    return ϵs_or_κs_new, exists_new, uniqs_new, approxs_new
end


function branch_continuation_helper_batch_fix_epsilon(
    ϵs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    exists = similar(uniqs)

    Threads.@threads for i in eachindex(ϵs, uniqs)
        ϵ = Arb(ϵs[i])
        Arblib.nonnegative_part!(ϵ, ϵ)

        G_x = x -> G(x..., ϵ, ξ₁, λ)
        dG_x = x -> G_jacobian_kappa(x..., ϵ, ξ₁, λ)

        exists[i] = verify_and_refine_root(G_x, dG_x, uniqs[i])
    end

    return exists
end

function branch_continuation_helper_batch_fix_kappa(
    κs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    exists = similar(uniqs)

    Threads.@threads for i in eachindex(κs, uniqs)
        κ = Arb(κs[i])

        G_x = x -> G(x[1:3]..., κ, x[4], ξ₁, λ)
        dG_x = x -> G_jacobian_kappa(x[1:3]..., κ, x[4], ξ₁, λ)

        exists[i] = verify_and_refine_root(G_x, dG_x, uniqs[i])
    end

    return exists
end

function branch_continuation_helper(
    ϵs_or_κs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    fix_kappa = false,
    pool = Distributed.WorkerPool(Distributed.workers()),
    batch_size = 32,
    verbose = false,
    log_progress = verose,
)
    # If the batch_size is too large then not all workers will get a
    # share. Take it smaller if this is the case.
    batch_size = min(batch_size, ceil(Int, length(ϵs_or_κs) / length(pool)))

    indices = firstindex(ϵs_or_κs):batch_size:lastindex(ϵs_or_κs)

    verbose && @info "Starting $(length(indices)) batch jobs of size $batch_size"

    tasks = map(indices) do index
        indices_batch = index:min(index + batch_size - 1, lastindex(ϵs_or_κs))

        if !fix_kappa
            @async Distributed.remotecall_fetch(
                (args...; kwargs...) -> @time(
                    branch_continuation_helper_batch_fix_epsilon(args...; kwargs...)
                ),
                pool,
                ϵs_or_κs[indices_batch],
                uniqs[indices_batch],
                ξ₁,
                λ,
            )
        else
            @async Distributed.remotecall_fetch(
                (args...; kwargs...) -> @time(
                    branch_continuation_helper_batch_fix_kappa(args...; kwargs...)
                ),
                pool,
                ϵs_or_κs[indices_batch],
                uniqs[indices_batch],
                ξ₁,
                λ,
            )
        end
    end

    verbose && @info "Collecting batch jobs"

    batches = fetch_with_progress(tasks, log_progress)

    exists = foldl(vcat, batches)

    return exists
end

function branch_continuation(
    ϵs_or_κs::Vector{NTuple{2,Arf}},
    exists::Vector{SVector{4,Arb}},
    uniqs::Vector{SVector{4,Arb}},
    approxs::Vector{SVector{4,Arb}},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    fix_kappa = false,
    pool = Distributed.WorkerPool(Distributed.workers()),
    batch_size = 32,
    maxevals::Integer = typemax(Int), # IMPROVE: Tune this
    depth::Integer = 20,
    verbose = false,
    log_progress = verbose,
)
    @assert length(ϵs_or_κs) == length(exists) == length(uniqs) == length(approxs)

    left_continuation_finite, left_continuation = check_left_continuation(exists, uniqs)

    failures_nonfinite = sum(!, left_continuation_finite)
    failures_containment = sum(!, left_continuation) - failures_nonfinite

    verbose && @info "Checking continuation for $(length(ϵs_or_κs)) subintervals"
    verbose && @info "Failures due to non-finite enclosures: $(failures_nonfinite)"
    verbose && @info "Failures due to containment: $(failures_containment)"

    if iszero(failures_containment)
        verbose && @info "No subintervals to bisect"

        return left_continuation, ϵs_or_κs, exists, uniqs, approxs
    end

    # We want to bisect all subintervals for which the enclosures are
    # finite but continuation fails. We want to bisect both the
    # subinterval itself, but also its left neighbour.
    to_bisect = left_continuation_finite .& .!left_continuation
    to_bisect[1:end-1] .|= to_bisect[2:end]

    verbose && @info "Bisecting $(sum(to_bisect)) subintervals"

    iterations = 0
    evals = 0

    verbose && @info "iteration: $(lpad(iterations, 2)), " *
          " starting intervals: $(lpad(2sum(to_bisect), 3))"

    while any(to_bisect) && !iszero(maxevals)
        iterations += 1
        evals += 2sum(to_bisect)

        ϵs_or_κs_bisected = ArbExtras.bisect_intervals(ϵs_or_κs, to_bisect)
        exists_bisected = permutedims(hcat(exists[to_bisect], exists[to_bisect]))[:]
        uniqs_bisected = permutedims(hcat(uniqs[to_bisect], uniqs[to_bisect]))[:]
        approxs_bisected = permutedims(hcat(approxs[to_bisect], approxs[to_bisect]))[:]

        exists_bisected_new = branch_continuation_helper(
            ϵs_or_κs_bisected,
            uniqs_bisected,
            ξ₁,
            λ;
            fix_kappa,
            pool,
            batch_size,
            verbose,
            log_progress,
        )

        if any(x -> !all(isfinite, x), exists_bisected_new)
            num_non_finite = count(x -> !all(isfinite, x), exists_bisected_new)
            verbose && @warn "Got $num_non_finite non-finite enclosures after bisection"

            # Reuse the old existence in this case. It is certainly
            # correct, but pessimistic.
            for i in findall(x -> !all(isfinite, x), exists_bisected_new)
                exists_bisected_new[i] = exists_bisected[i]
            end
        end

        # Insert the newly computed values back into the full vector
        ϵs_or_κs, exists, uniqs, approxs = continuation_insert_bisected(
            ϵs_or_κs,
            exists,
            uniqs,
            approxs,
            to_bisect,
            ϵs_or_κs_bisected,
            exists_bisected_new,
            uniqs_bisected,
            approxs_bisected,
        )

        left_continuation_finite, left_continuation = check_left_continuation(exists, uniqs)
        to_bisect = left_continuation_finite .& .!left_continuation
        to_bisect[1:end-1] .|= to_bisect[2:end]

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(2sum(to_bisect), 3))"

        if any(to_bisect) && evals >= maxevals
            verbose && @warn "Reached max evaluations" evals maxevals
            break
        end
        if any(to_bisect) && iterations >= depth
            verbose && @warn "Reached max depth" iterations depth
            break
        end
    end

    if any(to_bisect)
        verbose && @warn "Failed bisecting all subintervals for continuation"
    else
        verbose && @info "Succesfully bisected all subintervals for continuation!"
    end

    return left_continuation, ϵs_or_κs, exists, uniqs, approxs
end
