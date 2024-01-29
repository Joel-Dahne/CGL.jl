function verify_branch_segment_existence(
    (ϵ₁, ϵ₂)::NTuple{2,Arf},
    (μ₁, μ₂)::NTuple{2,Arb},
    (κ₁, κ₂)::NTuple{2,Arb},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    depth_start::Integer = ceil(Int, log2(Threads.nthreads())),
    maxevals::Integer = 1000,
    depth::Integer = 20,
    verbose = false,
)
    if ϵ₁ > ϵ₂
        error("TODO: Handle this case")
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

        Threads.@threads for i in eachindex(ϵs_remaining)
            ϵ = Arb(ϵs_remaining[i])
            Arblib.nonnegative_part!(ϵ, ϵ)
            λ_ϵ = CGLParams(λ; ϵ)

            # Linearly interpolate and refine
            μ, γ, κ = refine_approximation_with_interpolation(
                (μ₁, μ₂),
                (κ₁, κ₂),
                (Arb(ϵ₁), Arb(ϵ₂)),
                ξ₁,
                λ_ϵ,
            )

            # Minimum radius of κ for which we expect to be able to
            # prove existence
            r_min = let
                # Linearly interpolate κ at endpoints of ϵ.
                tₗ = (ϵ₂ - ϵs_remaining[i][1]) / (ϵ₂ - ϵ₁)
                tᵤ = (ϵ₂ - ϵs_remaining[i][2]) / (ϵ₂ - ϵ₁)

                κₗ_estimate = (1 - tₗ) * κ₁ + tₗ * κ₂
                κᵤ_estimate = (1 - tᵤ) * κ₁ + tᵤ * κ₂

                abs(κₗ_estimate - κᵤ_estimate)
            end

            # IMPROVE: Tune which values to try
            rs = [3, 1.8, 1.6, 1.4, 1.2, 1.1] .* r_min

            approxs_iteration[i] = SVector(μ, real(γ), imag(γ), κ)

            # Try solving
            exists_iteration[i], uniqs_iteration[i] = CGL.G_solve(
                μ,
                real(γ),
                imag(γ),
                κ,
                ξ₁,
                λ_ϵ,
                return_uniqueness = Val{true}(),
                verbose = false;
                rs,
            )
        end

        # Find all intervals for which existence was proved
        finished = map(x -> all(isfinite, x), exists_iteration)

        # Add all finished intervals to the result
        append!(ϵs_finished, ϵs_remaining[finished])
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
