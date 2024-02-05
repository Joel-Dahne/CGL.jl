function verify_branch_segment_existence_epsilon(
    (ϵ₁, ϵ₂)::NTuple{2,Arb},
    (μ₁, μ₂)::NTuple{2,Arb},
    (κ₁, κ₂)::NTuple{2,Arf},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    depth_start::Integer = ceil(Int, log2(Threads.nthreads())),
    maxevals::Integer = 1000,
    depth::Integer = 20,
    verbose = false,
)
    ArbExtras.check_interval(κ₂, κ₁)

    # List of intervals to check
    κs_remaining::Vector{NTuple{2,Arf}} =
        ArbExtras.bisect_interval_recursive(κ₁, κ₂, depth_start)

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

        Threads.@threads for i in eachindex(κs_remaining)
            κ = Arb(reverse(κs_remaining[i]))

            # Linearly interpolate and refine
            μ, γ, ϵ = refine_approximation_epsilon_with_interpolation(
                (μ₁, μ₂),
                (Arb(κ₁), Arb(κ₂)),
                (ϵ₁, ϵ₂),
                κ,
                ξ₁,
                λ,
            )

            # Minimum radius of ϵ for which we expect to be able to
            # prove existence
            r_min = let
                # Linearly interpolate ϵ at endpoints of κ.
                tₗ = (κ₂ - κs_remaining[i][1]) / (κ₂ - κ₁)
                tᵤ = (κ₂ - κs_remaining[i][2]) / (κ₂ - κ₁)

                ϵₗ_estimate = (1 - tₗ) * ϵ₁ + tₗ * ϵ₂
                ϵᵤ_estimate = (1 - tᵤ) * ϵ₁ + tᵤ * ϵ₂

                abs(ϵₗ_estimate - ϵᵤ_estimate)
            end

            # IMPROVE: Tune which values to try
            rs = [1.8, 1.4, 1.2] .* r_min

            approxs_iteration[i] = SVector(μ, real(γ), imag(γ), ϵ)

            # Try solving
            exists_iteration[i], uniqs_iteration[i] = CGL.G_solve_epsilon(
                μ,
                real(γ),
                imag(γ),
                κ,
                ξ₁,
                CGLParams(λ; ϵ),
                return_uniqueness = Val{true}(),
                verbose = false;
                rs,
            )
        end

        # Find all intervals for which existence was proved
        finished = map(x -> all(isfinite, x), exists_iteration)

        # Add all finished intervals to the result
        append!(κs_finished, κs_remaining[finished])
        append!(exists, exists_iteration[finished])
        append!(uniqs, uniqs_iteration[finished])
        append!(approxs, approxs_iteration[finished])

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
