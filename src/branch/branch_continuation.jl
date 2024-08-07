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
    N = length(ϵs_or_κs) + sum(to_bisect)

    ϵs_or_κs_new = typeof(ϵs_or_κs)(undef, N)
    exists_new = typeof(exists)(undef, N)
    uniqs_new = typeof(uniqs)(undef, N)
    approxs_new = typeof(approxs)(undef, N)

    j = 1
    k = 1
    for i in eachindex(ϵs_or_κs)
        if to_bisect[i]
            ϵs_or_κs_new[j:j+1] .= @view ϵs_or_κs_bisected[k:k+1]
            exists_new[j:j+1] .= @view exists_bisected[k:k+1]
            uniqs_new[j:j+1] .= @view uniqs_bisected[k:k+1]
            approxs_new[j:j+1] .= @view approxs_bisected[k:k+1]

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
    exists = tmap(eltype(uniqs), eachindex(ϵs, uniqs), scheduler = :greedy) do i
        ϵ = Arb(ϵs[i])
        Arblib.nonnegative_part!(ϵ, ϵ)

        G_x = x -> G(x..., ϵ, ξ₁, λ)
        dG_x = x -> G_jacobian_kappa(x..., ϵ, ξ₁, λ)

        verify_and_refine_root(G_x, dG_x, uniqs[i])
    end

    # There has been issues with high memory consumption giving
    # OOM crashes on SLURM. Explicitly galling gc here helps with
    # that.
    GC.gc()

    return exists
end

function branch_continuation_helper_batch_fix_kappa(
    κs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    exists = tmap(eltype(uniqs), eachindex(κs, uniqs), scheduler = :greedy) do i
        κ = Arb(κs[i])

        G_x = x -> G(x[1:3]..., κ, x[4], ξ₁, λ)
        dG_x = x -> G_jacobian_epsilon(x[1:3]..., κ, x[4], ξ₁, λ)

        verify_and_refine_root(G_x, dG_x, uniqs[i])
    end

    # There has been issues with high memory consumption giving
    # OOM crashes on SLURM. Explicitly galling gc here helps with
    # that.
    GC.gc()

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
    log_progress = verbose,
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

function branch_continuation_helper_G_solve_batch_fix_epsilon(
    ϵs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    exists::Vector{SVector{4,Arb}},
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    uniqs_new = similar(uniqs)
    exists_new = similar(exists)

    tforeach(eachindex(ϵs, uniqs, exists), scheduler = :greedy) do i
        ϵ = Arb(ϵs[i])
        Arblib.nonnegative_part!(ϵ, ϵ)

        # TODO: Try switching to the _alt version here. We have to be
        # a bit careful though. The main goal is not to expand the
        # uniqueness, but is it still worth it to set
        # try_expand_uniqueness to true? We might also want to refine
        # the approximation first? The midpoint for exists[i] is the
        # midpoint for the non-bisected interval and is not the best
        # approximation for the current interval.

        # It is fine to take a very small r since we only update the
        # uniqueness if it is actually larger than before. We are
        # mostly interested in improved enclosure for existence.
        rs = [1.2, 1.1, 1, 0.9, 0.8, 0.1, 1e-2] * radius(Arb, uniqs[i][4])

        exists_new[i], uniqs_new[i] = G_solve_fix_epsilon(
            midpoint.(Arb, exists[i])...,
            ϵ,
            ξ₁,
            λ,
            return_uniqueness = Val{true}();
            rs,
        )
    end

    # There has been issues with high memory consumption giving
    # OOM crashes on SLURM. Explicitly galling gc here helps with
    # that.
    GC.gc()

    return uniqs_new, exists_new
end

function branch_continuation_helper_G_solve_batch_fix_kappa(
    κs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    exists::Vector{SVector{4,Arb}},
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    uniqs_new = similar(uniqs)
    exists_new = similar(exists)

    tforeach(eachindex(κs, uniqs, exists), scheduler = :greedy) do i
        κ = Arb(κs[i])

        # TODO: Try switching to the _alt version here. We have to be
        # a bit careful though, as discussed in the _fix_epsilon
        # version above.

        # It is fine to take a very small r since we only update the
        # uniqueness if it is actually larger than before. We are
        # mostly interested in improved enclosure for existence.
        rs = [1.2, 1.1, 1, 0.9, 0.8, 0.1, 1e-2] * radius(Arb, uniqs[i][4])

        exists_new[i], uniqs_new[i] = G_solve_fix_kappa(
            midpoint.(Arb, exists[i][1:3])...,
            κ,
            midpoint(Arb, exists[i][4]),
            ξ₁,
            λ,
            return_uniqueness = Val{true}();
            rs,
        )
    end

    # There has been issues with high memory consumption giving
    # OOM crashes on SLURM. Explicitly galling gc here helps with
    # that.
    GC.gc()

    return uniqs_new, exists_new
end

function branch_continuation_G_solve_helper(
    ϵs_or_κs::Vector{NTuple{2,Arf}},
    uniqs::Vector{SVector{4,Arb}},
    exists::Vector{SVector{4,Arb}},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    fix_kappa = false,
    pool = Distributed.WorkerPool(Distributed.workers()),
    batch_size = 32,
    verbose = false,
    log_progress = verbose,
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
                    branch_continuation_helper_G_solve_batch_fix_epsilon(
                        args...;
                        kwargs...,
                    )
                ),
                pool,
                ϵs_or_κs[indices_batch],
                uniqs[indices_batch],
                exists[indices_batch],
                ξ₁,
                λ,
            )
        else
            @async Distributed.remotecall_fetch(
                (args...; kwargs...) -> @time(
                    branch_continuation_helper_G_solve_batch_fix_kappa(args...; kwargs...)
                ),
                pool,
                ϵs_or_κs[indices_batch],
                uniqs[indices_batch],
                exists[indices_batch],
                ξ₁,
                λ,
            )
        end
    end

    verbose && @info "Collecting batch jobs"

    batches = fetch_with_progress(tasks, log_progress)

    uniqs_new = foldl(vcat, getindex.(batches, 1))
    exists_new = foldl(vcat, getindex.(batches, 2))

    return uniqs_new, exists_new
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
    batch_size = 128,
    maxevals::Integer = typemax(Int), # IMPROVE: Tune this
    depth::Integer = 20,
    verbose = false,
    log_progress = verbose,
)
    @assert length(ϵs_or_κs) == length(exists) == length(uniqs) == length(approxs)

    # For the turn and the bottom part of the branch the intervals are
    # reversely sorted. For ArbExtras.bisect_intervals to work as
    # expected they need to be reversed in that case. We therefore
    # need to keep track of this.
    rev = issorted(ϵs_or_κs, by = interval -> interval[1], rev = true)

    rev ||
        issorted(ϵs_or_κs, by = interval -> interval[1]) ||
        throw(ArgumentError("expected input to be sorted"))

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
    # finite but continuation fails.
    # TODO: In the end we might want to only allow cases where all are
    # finite? Simplifies the code a bit.
    to_bisect = left_continuation_finite .& .!left_continuation

    verbose && @info "Bisecting $(sum(to_bisect)) subintervals"

    iterations = 0
    evals = 0

    verbose && @info "iteration: $(lpad(iterations, 2)), " *
          " starting intervals: $(lpad(2sum(to_bisect), 3))"

    while any(to_bisect) && !iszero(maxevals)
        iterations += 1
        evals += 2sum(to_bisect)

        if !rev
            ϵs_or_κs_bisected = ArbExtras.bisect_intervals(ϵs_or_κs, to_bisect)
        else
            # ArbExtras.bisect_intervals needs the input to be sorted
            # to give the expected results. In this case we therefore
            # reverse it twice
            reverse!(ϵs_or_κs)
            reverse!(to_bisect)
            ϵs_or_κs_bisected = reverse!(ArbExtras.bisect_intervals(ϵs_or_κs, to_bisect))
            reverse!(ϵs_or_κs)
            reverse!(to_bisect)
        end

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
            verbose &&
                @info "$num_non_finite subintervals failed verification after bisection"
            verbose && @info "Trying to verify them using G_solve"

            indices_failed = findall(x -> !all(isfinite, x), exists_bisected_new)

            uniqs_bisected_G_solve, exists_bisected_G_solve =
                branch_continuation_G_solve_helper(
                    ϵs_or_κs_bisected[indices_failed],
                    uniqs_bisected[indices_failed],
                    exists_bisected[indices_failed],
                    ξ₁,
                    λ;
                    fix_kappa,
                    pool,
                    batch_size,
                    verbose,
                    log_progress,
                )

            if verbose && any(x -> !all(isfinite, x), exists_bisected_G_solve)
                num_non_finite_G_solve =
                    count(x -> !all(isfinite, x), exists_bisected_G_solve)
                @warn "$num_non_finite_G_solve subintervals failed verification after G_solve"
            end

            for i in eachindex(indices_failed)
                if all(isfinite, exists_bisected_G_solve[i])
                    # Update the existence
                    exists_bisected_new[indices_failed[i]] = exists_bisected_G_solve[i]
                    # If the new interval for uniqueness is larger
                    # (which is typically the case) update that as
                    # well.
                    if all(
                        Arblib.contains.(
                            uniqs_bisected_G_solve[i],
                            uniqs_bisected[indices_failed[i]],
                        ),
                    )
                        uniqs_bisected[indices_failed[i]] = uniqs_bisected_G_solve[i]
                    end
                else
                    # Reuse the old existence in this case. It is
                    # certainly correct, but pessimistic.
                    exists_bisected_new[indices_failed[i]] =
                        exists_bisected[indices_failed[i]]
                end
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
