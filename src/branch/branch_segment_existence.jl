function branch_segment_existence_fix_epsilon(
    (μ₁, μ₂)::NTuple{2,Arb},
    (κ₁, κ₂)::NTuple{2,Arb},
    (ϵ₁, ϵ₂)::NTuple{2,Arf},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    depth_start::Integer = 5,
    maxevals::Integer = 1000,
    depth::Integer = 20,
    try_expand_uniqueness = true,
    verbose = false,
)
    if ϵ₁ > ϵ₂
        # Flip order of arguments and reverse result
        res = branch_segment_existence_fix_epsilon(
            (μ₂, μ₁),
            (κ₂, κ₁),
            (ϵ₂, ϵ₁),
            ξ₁,
            λ;
            depth_start,
            maxevals,
            depth,
            verbose,
        )

        return reverse!.(res)
    end

    ArbExtras.check_interval(ϵ₁, ϵ₂)

    # List of intervals to check
    ϵs_remaining::Vector{NTuple{2,Arf}} =
        ArbExtras.bisect_interval_recursive(ϵ₁, ϵ₂, depth_start)

    # Stores the result
    ϵs_finished = empty(ϵs_remaining)
    exists = empty(ϵs_remaining, SVector{4,Arb})
    uniqs = empty(ϵs_remaining, SVector{4,Arb})
    approxs = empty(ϵs_remaining, SVector{4,Arb})

    iterations = 0
    evals = 0

    verbose && @info "iteration: $(lpad(iterations, 2)), " *
          " starting intervals: $(lpad(length(ϵs_remaining), 3))"

    while !isempty(ϵs_remaining)
        iterations += 1
        evals += length(ϵs_remaining)

        exists_iteration = similar(ϵs_remaining, eltype(exists))
        uniqs_iteration = similar(ϵs_remaining, eltype(uniqs))
        approxs_iteration = similar(ϵs_remaining, eltype(approxs))

        tforeach(eachindex(ϵs_remaining), scheduler = :greedy) do i
            ϵ = Arb(ϵs_remaining[i])
            Arblib.nonnegative_part!(ϵ, ϵ)

            # Linearly interpolate and refine
            μ, γ, κ = refine_approximation_fix_epsilon_with_interpolation(
                (μ₁, μ₂),
                (κ₁, κ₂),
                (Arb(ϵ₁), Arb(ϵ₂)),
                ϵ,
                ξ₁,
                λ,
            )

            approxs_iteration[i] = SVector(μ, real(γ), imag(γ), κ)

            exists_iteration[i], uniqs_iteration[i] = CGL.G_solve_fix_epsilon(
                μ,
                real(γ),
                imag(γ),
                κ,
                ϵ,
                ξ₁,
                λ,
                return_uniqueness = Val{true}(),
                expansion_rate = 0.1;
                try_expand_uniqueness,
            )
        end

        # Find all intervals for which existence was proved
        finished = map(x -> all(isfinite, x), exists_iteration)

        # Add all finished intervals to the result
        append!(ϵs_finished, ϵs_remaining[finished])
        append!(exists, exists_iteration[finished])
        append!(uniqs, uniqs_iteration[finished])
        append!(approxs, approxs_iteration[finished])

        # There has been issues with high memory consumption giving
        # OOM crashes on SLURM. Explicitly galling gc here helps with
        # that.
        GC.gc()

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(sum(.!finished), 3))"

        # Check if we have done the maximum number of evaluations or
        # reached the maximum depth
        if !all(finished) && (evals >= maxevals || iterations >= depth - depth_start)
            if verbose
                evals >= maxevals &&
                    @info "reached maximum number of evaluations $evals >= $maxevals"
                iterations >= depth - depth_start && @info "reached maximum depth $depth"
            end

            # Append also the unfinished intervals to the result
            append!(ϵs_finished, ϵs_remaining[.!finished])
            append!(exists, exists_iteration[.!finished])
            append!(uniqs, uniqs_iteration[.!finished])
            append!(approxs, approxs_iteration[.!finished])

            # Don't bisect any more intervals
            ϵs_remaining = empty(ϵs_remaining)
        else
            # Bisect the remaining intervals
            ϵs_remaining = ArbExtras.bisect_intervals(ϵs_remaining, .!finished)
        end
    end

    p = sortperm(ϵs_finished, by = interval -> interval[1])
    return ϵs_finished[p], exists[p], uniqs[p], approxs[p]
end

function branch_segment_existence_fix_kappa(
    (μ₁, μ₂)::NTuple{2,Arb},
    (κ₁, κ₂)::NTuple{2,Arf},
    (ϵ₁, ϵ₂)::NTuple{2,Arb},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    depth_start::Integer = 5,
    maxevals::Integer = 1000,
    depth::Integer = 20,
    try_expand_uniqueness = true,
    verbose = false,
)
    # We should have κ₂ < κ₁
    ArbExtras.check_interval(κ₂, κ₁)

    # List of intervals to check
    κs_remaining::Vector{NTuple{2,Arf}} =
        ArbExtras.bisect_interval_recursive(κ₂, κ₁, depth_start)

    # Stores the result
    κs_finished = empty(κs_remaining)
    exists = empty(κs_remaining, SVector{4,Arb})
    uniqs = empty(κs_remaining, SVector{4,Arb})
    approxs = empty(κs_remaining, SVector{4,Arb})

    iterations = 0
    evals = 0

    verbose && @info "iteration: $(lpad(iterations, 2)), " *
          " starting intervals: $(lpad(length(κs_remaining), 3))"

    while !isempty(κs_remaining)
        iterations += 1
        evals += length(κs_remaining)

        exists_iteration = similar(κs_remaining, eltype(exists))
        uniqs_iteration = similar(κs_remaining, eltype(uniqs))
        approxs_iteration = similar(κs_remaining, eltype(approxs))

        tforeach(eachindex(κs_remaining), scheduler = :greedy) do i
            κ = Arb(κs_remaining[i])

            # Linearly interpolate and refine
            μ, γ, ϵ = refine_approximation_fix_kappa_with_interpolation(
                (μ₁, μ₂),
                (Arb(κ₁), Arb(κ₂)),
                (ϵ₁, ϵ₂),
                κ,
                ξ₁,
                λ,
            )

            approxs_iteration[i] = SVector(μ, real(γ), imag(γ), ϵ)

            exists_iteration[i], uniqs_iteration[i] = CGL.G_solve_fix_kappa(
                μ,
                real(γ),
                imag(γ),
                κ,
                ϵ,
                ξ₁,
                λ,
                return_uniqueness = Val{true}(),
                expansion_rate = 0.1;
                try_expand_uniqueness,
            )
        end

        # Find all intervals for which existence was proved
        finished = map(x -> all(isfinite, x), exists_iteration)

        # Add all finished intervals to the result
        append!(κs_finished, κs_remaining[finished])
        append!(exists, exists_iteration[finished])
        append!(uniqs, uniqs_iteration[finished])
        append!(approxs, approxs_iteration[finished])

        # There has been issues with high memory consumption giving
        # OOM crashes on SLURM. Explicitly galling gc here helps with
        # that.
        GC.gc()

        verbose && @info "iteration: $(lpad(iterations, 2)), " *
              "remaining intervals: $(lpad(sum(.!finished), 3))"

        # Check if we have done the maximum number of evaluations or
        # reached the maximum depth
        if !all(finished) && (evals >= maxevals || iterations >= depth - depth_start)
            if verbose
                evals >= maxevals &&
                    @info "reached maximum number of evaluations $evals >= $maxevals"
                iterations >= depth - depth_start && @info "reached maximum depth $depth"
            end

            # Append also the unfinished intervals to the result
            append!(κs_finished, κs_remaining[.!finished])
            append!(exists, exists_iteration[.!finished])
            append!(uniqs, uniqs_iteration[.!finished])
            append!(approxs, approxs_iteration[.!finished])

            # Don't bisect any more intervals
            κs_remaining = empty(κs_remaining)
        else
            # Bisect the remaining intervals
            κs_remaining = ArbExtras.bisect_intervals(κs_remaining, .!finished)
        end
    end

    p = sortperm(κs_finished, by = interval -> interval[1], rev = true)
    return κs_finished[p], exists[p], uniqs[p], approxs[p]
end
