function _bisect_and_fill(
    ϵs::Vector{NTuple{2,Arb}},
    exists::Vector{SVector{4,Arb}},
    uniqs::Vector{SVector{4,Arb}},
    to_bisect::Vector{Bool},
)
    N = length(ϵs) + sum(to_bisect) # New length
    ϵs_new = Vector{eltype(ϵs)}(undef, N)
    exists_new = Vector{eltype(exists)}(undef, N)
    uniqs_new = Vector{eltype(uniqs)}(undef, N)
    to_bisect_new = Vector{eltype(to_bisect)}(undef, N)

    indeterminate_vector = SVector(
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
    )

    j = 1
    for i in eachindex(ϵs, exists, uniqs, to_bisect)
        if to_bisect[i]
            a, b = ϵs[i]
            mid = ArbExtras._midpoint_interval(a, b)

            ϵs_new[j], ϵs_new[j+1] = (a, mid), (mid, b)
            # Existence and uniqueness propagates to the subintervals
            exists_new[j], exists_new[j+1] = copy(exists[i]), copy(exists[i])
            uniqs_new[j], uniqs_new[j+1] = copy(uniqs[i]), copy(uniqs[i])
            to_bisect_new[j], to_bisect_new[j+1] = true, true

            j += 2
        else
            ϵs_new[j] = ϵs[i]
            exists_new[j] = exists[i]
            uniqs_new[j] = uniqs[i]
            to_bisect_new[j] = to_bisect[i]

            j += 1
        end
    end

    return ϵs_new, exists_new, uniqs_new, to_bisect_new
end

function check_continuation(exists::Vector{SVector{4,Arb}}, uniqs::Vector{SVector{4,Arb}})
    to_bisect = fill(false, length(exists))

    for i in Iterators.drop(eachindex(to_bisect, exists, uniqs), 1)
        continuation_ok =
            all(isfinite, exists[i-1]) &&
            all(isfinite, uniqs[i]) &&
            all(Arblib.contains.(uniqs[i-1], exists[i]))
        if !continuation_ok
            to_bisect[i-1] = true
            to_bisect[i] = true
        end
    end

    return to_bisect
end

function verify_branch_segment(
    (ϵ₁, ϵ₂)::Tuple{Arb,Arb},
    (μ₁, μ₂)::Tuple{Arb,Arb},
    (κ₁, κ₂)::Tuple{Arb,Arb},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    verbose = false,
)
    reversed = ϵ₁ > ϵ₂
    if reversed
        ϵ₁, ϵ₂ = ϵ₂, ϵ₁
        μ₁, μ₂ = μ₂, μ₁
        κ₁, κ₂ = κ₂, κ₁
    end

    # Refine approximations at endpoints and find corresponding γ
    λ₁ = CGLParams(λ, ϵ = ϵ₁)
    λ₂ = CGLParams(λ, ϵ = ϵ₂)

    μ₁, γ₁, κ₁ = refine_approximation(μ₁, κ₁, ξ₁, λ₁)
    μ₂, γ₂, κ₂ = refine_approximation(μ₂, κ₂, ξ₁, λ₂)

    # Used for representing enclosures that are yet to be computed
    indeterminate_vector = SVector(
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
    )

    # Vector of subintervals of (ϵ₁, ϵ₂), represented by tuples of Arb
    ϵs = [(ϵ₁, ϵ₂)]
    # For each subinterval these holds boxes of existence and
    # uniqueness of the solutions.
    exists = SVector{4,Arb}[indeterminate_vector]
    uniqs = SVector{4,Arb}[indeterminate_vector]
    # Holds information about which subintervals should be bisected
    # for the next iteration
    to_bisect = [true]

    # Bisect a few times to make use of all threads from the start
    while 2length(ϵs) <= Threads.nthreads()
        ϵs, exists, uniqs, to_bisect = _bisect_and_fill(ϵs, exists, uniqs, to_bisect)
    end

    if verbose
        @info "Bisecting for finite enclosures"
        @info "iteration: $(lpad(0, 2)), starting intervals: $(lpad(sum(to_bisect), 3))"
    end

    max_iterations = 20
    max_evals = 1000
    iteration = 0
    evals = 0

    while any(to_bisect)
        iteration += 1
        evals += sum(to_bisect)

        last_iteration = iteration >= max_iterations || evals >= max_evals

        Threads.@threads for i in findall(to_bisect)
            ϵ = Arb(ϵs[i])
            Arblib.nonnegative_part!(ϵ, ϵ)
            λ_ϵ = CGLParams(λ; ϵ)

            # Linearly interpolate and refine
            μ, γ, κ = refine_approximation_with_interpolation(
                (μ₁, μ₂),
                (κ₁, κ₂),
                (ϵ₁, ϵ₂),
                ξ₁,
                λ_ϵ,
            )

            # Minimum radius of κ for which we expect to be able to
            # prove existence
            r_min = let
                # Linearly interpolate κ at endpoints of ϵ.
                ϵₗ, ϵᵤ = ϵs[i]
                tₗ = (ϵ₂ - ϵₗ) / (ϵ₂ - ϵ₁)
                tᵤ = (ϵ₂ - ϵᵤ) / (ϵ₂ - ϵ₁)

                κₗ_estimate = (1 - tₗ) * κ₁ + tₗ * κ₂
                κᵤ_estimate = (1 - tᵤ) * κ₁ + tᵤ * κ₂

                abs(κₗ_estimate - κᵤ_estimate)
            end

            # IMPROVE: Tune which values to try
            rs = [3, 1.8, 1.6, 1.4, 1.2, 1.1] .* r_min

            # Try solving
            exist, uniq = CGL.G_solve(
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

            if last_iteration || all(isfinite, exist)
                # Add to result and don't split further
                exists[i] = exist
                uniqs[i] = uniq
                to_bisect[i] = false
            else
                # Split further
                to_bisect[i] = true
            end
        end

        ϵs, exists, uniqs, to_bisect = _bisect_and_fill(ϵs, exists, uniqs, to_bisect)

        verbose && @info "iteration: $(lpad(iteration, 2)), " *
              "remaining intervals: $(lpad(sum(to_bisect) ÷ 2, 3))"
    end

    if verbose
        failures_existence = count(x -> !all(isfinite, x), exists)

        failures_existence > 0 &&
            @warn "Could not prove existence on all subintervals" failures_existence
    end

    # Continue to bisect intervals where unique continuation is not determined
    to_bisect = check_continuation(exists, uniqs)

    if verbose && any(to_bisect)
        @info "Bisecting for continuation"
        @info "iteration: $(lpad(0, 2)), remaining intervals: $(lpad(sum(to_bisect), 3))"
    end

    ϵs, exists, uniqs, to_bisect = _bisect_and_fill(ϵs, exists, uniqs, to_bisect)

    max_iterations = 20
    max_evals = 1000
    iteration = 0
    evals = 0
    last_iteration = false

    while !last_iteration && any(to_bisect)
        iteration += 1
        evals += sum(to_bisect)

        last_iteration = iteration >= max_iterations || evals >= max_evals
        failures = 0

        Threads.@threads for i in findall(to_bisect)
            ϵ = Arb(ϵs[i])
            Arblib.nonnegative_part!(ϵ, ϵ)
            λ_ϵ = CGLParams(λ; ϵ)

            # Linearly interpolate and refine
            μ, γ, κ = refine_approximation_with_interpolation(
                (μ₁, μ₂),
                (κ₁, κ₂),
                (ϵ₁, ϵ₂),
                ξ₁,
                λ_ϵ,
            )

            # Minimum radius of κ for which we expect to be able to
            # prove continuation
            r_min = let
                # Linearly interpolate κ at endpoints one subinterval
                # away, assuming it has the same radius.
                tₗ = t = (ϵ₂ - (ϵₗ - 2radius(ϵ))) / (ϵ₂ - ϵ₁)
                tᵤ = t = (ϵ₂ - (ϵᵤ + 2radius(ϵ))) / (ϵ₂ - ϵ₁)

                κₗ_estimate = (1 - tₗ) * κ₁ + tₗ * κ₂
                κᵤ_estimate = (1 - tᵤ) * κ₁ + tᵤ * κ₂

                abs(κₗ_estimate - κᵤ_estimate)
            end

            # IMPROVE: Tune which values to try
            rs = [3, 1.8, 1.6, 1.4, 1.2, 1.1] .* r_min

            # Try solving
            exist, uniq = CGL.G_solve(
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

            if !all(isfinite, exist) && all(isfinite, exists[i])
                # We got a finite result earlier, but not anymore
                failures += 1
            else
                # Add to results
                exists[i] = exist
                uniqs[i] = uniq
            end
        end

        to_bisect = check_continuation(exists, uniqs)

        verbose && @info "iteration: $(lpad(iteration, 2)), " *
              "remaining intervals: $(lpad(sum(to_bisect), 3))" *
              ifelse(iszero(failures), "", ", failures: $(lpad(failures, 3))")

        if !last_iteration
            ϵs, exists, uniqs, to_bisect = _bisect_and_fill(ϵs, exists, uniqs, to_bisect)
        end
    end

    to_bisect = check_continuation(exists, uniqs)

    verbose &&
        any(to_bisect) &&
        @warn "Could not prove continuation on all subintervals" sum(to_bisect)

    if reversed
        reverse!(ϵs)
        reverse!(exists)
        reverse!(uniqs)
    end

    return ϵs, exists, uniqs
end

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

function verify_branch(
    ϵs::Vector{Arb},
    μs::Vector{Arb},
    κs::Vector{Arb},
    ξ₁::Arb,
    λ::CGLParams{Arb};
    verbose = false,
    verbose_segments = false,
    log_progress = false,
)
    @assert length(ϵs) == length(μs) == length(κs)

    pool = Distributed.WorkerPool(Distributed.workers())

    verbose && @info "Verifying $(length(ϵs)-1) segments"

    tasks = map(1:length(ϵs)-1) do i
        @async Distributed.remotecall_fetch(
            verify_branch_segment,
            pool,
            (ϵs[i], ϵs[i+1]),
            (μs[i], μs[i+1]),
            (κs[i], κs[i+1]),
            ξ₁,
            λ,
            verbose = verbose_segments,
        )
    end

    if log_progress
        @progress segments = [fetch(task) for task in tasks]
    else
        segments = [fetch(task) for task in tasks]
    end

    failed_segments = map(segments) do (_, exists, uniqs)
        any(check_continuation(exists, uniqs))
    end

    if iszero(any(failed_segments))
        verbose && @info "Succesfully verified all segments"
    else
        verbose && @warn "Failed verifying $(sum(failed_segments)) segments"
    end

    ϵs = reduce(vcat, getindex.(segments, 1))
    exists = reduce(vcat, getindex.(segments, 2))
    uniqs = reduce(vcat, getindex.(segments, 3))

    to_bisect = check_continuation(exists, uniqs)

    if iszero(any(failed_segments))
        # Check if any endpoints needs to be bisected
        if any(to_bisect)
            verbose &&
                @info "Continuation failed for $(sum(to_bisect)) intersections of segments"

            # TODO: Implement this
            verbose && @error "Bisection for intersections not yet implemented"
        else
            verbose && @info "Verified continuation for all intersections of segments"
        end
    else
        verbose && @warn "Not attempting to veryify intersection of segments"
    end

    success = .!check_continuation(exists, uniqs)

    return success, ϵs, exists, uniqs
end
