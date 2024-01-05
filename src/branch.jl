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

    fill_value = SVector(
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
            exists_new[j], exists_new[j+1] = fill_value, fill_value
            uniqs_new[j], uniqs_new[j+1] = fill_value, fill_value
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
    (ϵ₁, ϵ₂),
    (μ₁, μ₂),
    (κ₁, κ₂),
    ξ₁,
    λ::CGLParams;
    verbose = false,
)
    reversed = ϵ₁ > ϵ₂
    if reversed
        ϵ₁, ϵ₂ = ϵ₂, ϵ₁
        μ₁, μ₂ = μ₂, μ₁
        κ₁, κ₂ = κ₂, κ₁
    end

    # Refine approximations and endpoints and find corresponding γ
    λ₁ = CGLParams(λ, ϵ = ϵ₁)
    λ₂ = CGLParams(λ, ϵ = ϵ₂)

    μ₁, γ₁, κ₁ = refine_approximation(μ₁, κ₁, ξ₁, λ₁)
    μ₂, γ₂, κ₂ = refine_approximation(μ₂, κ₂, ξ₁, λ₂)

    if false
        # Check that the endpoints can be verified
        r1 = CGL.G_solve(μ₁, γ₁, κ₁, ξ₁, λ₁)
        r2 = CGL.G_solve(μ₂, γ₂, κ₂, ξ₁, λ₂)

        @assert all(isfinite, r1)
        @assert all(isfinite, r2)
    end

    # Used for representing enclosures that are yet to be computed
    fill_value = SVector(
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
    )

    # Vector of subintervals of (ϵ₁, ϵ₂), represented by tuples of Arb
    ϵs = [(ϵ₁, ϵ₂)]
    # For each subinterval these holds boxes of existence and
    # uniqueness of the solutions.
    exists = SVector{4,Arb}[fill_value]
    uniqs = SVector{4,Arb}[fill_value]
    # Holds information about which subintervals should be bisected
    # for the next iteration
    to_bisect = [true]

    # Bisect a few times to make use of all threads from the start
    while length(ϵs) < Threads.nthreads()
        ϵs, exists, uniqs, to_bisect = _bisect_and_fill(ϵs, exists, uniqs, to_bisect)
    end

    verbose && @info "Bisecting for finite enclosures"

    max_iterations = 10
    iteration = 0

    verbose && @info "iteration: $(lpad(iteration, 2)), " *
          "remaining intervals: $(lpad(sum(to_bisect), 3))"

    while any(to_bisect)
        iteration += 1

        Threads.@threads for i in findall(to_bisect)
            ϵ = Arb(ϵs[i])
            Arblib.nonnegative_part!(ϵ, ϵ)
            λ_ϵ = CGLParams(λ; ϵ)

            # Linearly interpolate as a first approximation
            t = (ϵ₂ - midpoint(ϵ)) / (ϵ₂ - ϵ₁)

            μ = (1 - t) * μ₁ + t * μ₂
            γ = (1 - t) * γ₁ + t * γ₂
            κ = (1 - t) * κ₁ + t * κ₂

            # Refine approximation
            μ, γ, κ = refine_approximation(μ, γ, κ, ξ₁, λ_ϵ)

            # Minimum r for which we expect to be able to prove
            # continuation of the branch
            r_min = let
                ϵₗ, ϵᵤ = ϵs[i]
                # Estimated κ values at ϵₗ - (ϵᵤ - ϵₗ) and ϵᵤ + (ϵᵤ -
                # ϵₗ) using linear interpolation between ϵ₁ and ϵ₂
                tₗ = t = (ϵ₂ - (ϵₗ - (ϵᵤ - ϵₗ))) / (ϵ₂ - ϵ₁)
                tᵤ = t = (ϵ₂ - (ϵᵤ + (ϵᵤ - ϵₗ))) / (ϵ₂ - ϵ₁)
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

            if all(isfinite, exist)
                # Add to result and don't split further
                exists[i] = exist
                uniqs[i] = uniq
                to_bisect[i] = false
            elseif iteration >= max_iterations
                error("What to do?")
            else
                # Split further
                to_bisect[i] = true
            end
        end

        ϵs, exists, uniqs, to_bisect = _bisect_and_fill(ϵs, exists, uniqs, to_bisect)

        verbose && @info "iteration: $(lpad(iteration, 2)), " *
              "remaining intervals: $(lpad(sum(to_bisect), 3))"
    end

    verbose && @info "Bisecting for continuation"

    # Continue to bisect intervals where unique continuation is not determined
    to_bisect = check_continuation(exists, uniqs)

    iteration = 0

    verbose && @info "iteration: $(lpad(iteration, 2)), " *
          "remaining intervals: $(lpad(sum(to_bisect), 3))"

    while any(to_bisect)
        iteration += 1

        ϵs, exists, uniqs, to_bisect = _bisect_and_fill(ϵs, exists, uniqs, to_bisect)

        failures = 0

        Threads.@threads for i in findall(to_bisect)
            ϵ = Arb(ϵs[i])
            Arblib.nonnegative_part!(ϵ, ϵ)
            λ_ϵ = CGLParams(λ; ϵ)

            # Linearly interpolate as a first approximation
            t = (ϵ₂ - midpoint(ϵ)) / (ϵ₂ - ϵ₁)

            μ = (1 - t) * μ₁ + t * μ₂
            γ = (1 - t) * γ₁ + t * γ₂
            κ = (1 - t) * κ₁ + t * κ₂

            # Refine approximation
            μ, γ, κ = refine_approximation(μ, γ, κ, ξ₁, λ_ϵ)

            # Minimum r for which we expect to be able to prove
            # continuation of the branch
            r_min = let
                ϵₗ, ϵᵤ = ϵs[i]
                # Estimated κ values at ϵₗ - (ϵᵤ - ϵₗ) and ϵᵤ + (ϵᵤ -
                # ϵₗ) using linear interpolation between ϵ₁ and ϵ₂
                tₗ = t = (ϵ₂ - (ϵₗ - (ϵᵤ - ϵₗ))) / (ϵ₂ - ϵ₁)
                tᵤ = t = (ϵ₂ - (ϵᵤ + (ϵᵤ - ϵₗ))) / (ϵ₂ - ϵ₁)
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

            if !all(isfinite, exist) && all(isfinite, uniq)
                failures += 1
            end

            # Add to results
            exists[i] = exist
            uniqs[i] = uniq
        end

        to_bisect = check_continuation(exists, uniqs)

        verbose && @info "iteration: $(lpad(iteration, 2)), " *
              "remaining intervals: $(lpad(sum(to_bisect), 3))" *
              ifelse(iszero(failures), "", ", failures: $(lpad(failures, 3))")
    end

    if reversed
        reverse!(ϵs)
        reverse!(exists)
        reverse!(uniqs)
    end

    return ϵs, exists, uniqs
end

function verify_branch_points(ϵs, μs, κs, ξ₁, λ::CGLParams; verbose = false)
    res = similar(ϵs, SVector{4,Arb})

    p = ProgressMeter.Progress(length(res), enabled = verbose)

    Threads.@threads for i in eachindex(ϵs, μs, κs)
        λ_ϵ = CGLParams(λ, ϵ = ϵs[i])

        # Refine approximation
        μ, γ, κ = refine_approximation(μs[i], κs[i], ξ₁, λ_ϵ)

        res[i] = CGL.G_solve(μ, real(γ), imag(γ), κ, ξ₁, λ_ϵ)

        ProgressMeter.next!(p)
    end

    ProgressMeter.finish!(p)

    return res
end

function verify_branch(ϵs, μs, κs, ξ₁, λ::CGLParams; verbose = false)
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
            λ;
            verbose,
        )
    end

    segments = ProgressMeter.@showprogress map(fetch, tasks)

    failures = sum(segments) do (ϵs, exists, uniqs)
        any(check_continuation(exists, uniqs))
    end

    if iszero(failures)
        verbose && @info "Succesfully verified all segments"
    else
        verbose && @warn "Failed verifying $failues segments of $(length(segments))"
    end

    ϵs = reduce(vcat, getindex.(segments, 1))
    exists = reduce(vcat, getindex.(segments, 2))
    uniqs = reduce(vcat, getindex.(segments, 3))

    # Check if any endpoints needs to be bisected
    to_bisect = check_continuation(exists, uniqs)

    verbose && @info "Verifying intersections of segments" sum(to_bisect)

    return ϵs, exists, uniqs
end
