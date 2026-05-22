"""
count_critical_points(μ::Arb, γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}; verbose)
"""
function count_critical_points(
    μ::Arb,
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb};
    verbose = false,
)
    # Verify monotonicity on (ξ₁, ∞)
    monotone_ξ₁_inf = verify_monotonicity_infinity(γ, κ, ϵ, ξ₁, Λ; verbose)

    if monotone_ξ₁_inf
        verbose && @info "Verified monotonicity on (ξ₁, ∞)"
    else
        verbose && @warn "Could not verify monotonicity on (ξ₁, ∞)"
        return false, Arb[], Bool[]
    end

    # Compute enclosure of the profile on [0, ξ₁]

    ξs, Qs, d2Qs, abs2_Q_derivatives, abs2_Q_derivative2s =
        Q_zero_capd_curve(μ, κ, ϵ, ξ₁, Λ)

    if !all(Q -> all(isfinite, Q), Qs)
        verbose && @warn "Could not enclose curve on [0, ξ₁]"
        return false, Arb[], Bool[]
    end

    # Pick ξ₀ value
    i = findfirst(!Arblib.contains_zero, abs2_Q_derivatives)
    ξ₀ = ubound(Arb, ξs[i-1])
    @assert i >= 2

    # Verify monotonicity on [0, ξ₀]
    verified_zero = if Λ.d != 3
        # Just check that second order derivative is non-zero
        all(!Arblib.contains_zero, abs2_Q_derivative2s[1:(i-1)])
    else
        # For d = 3 we get very bad bounds for d2Qs near zero. Instead
        # we evaluate the Taylor expansion directly on [0, ξ₀] to
        # bound it
        (a_ξ₀, b_ξ₀, α_ξ₀, β_ξ₀), (d2a_ξ₀, d2b_ξ₀) =
            Q_zero_taylor(μ, κ, ϵ, ξ₀, Λ, enclose_curve = Val{true}())

        abs2_Q_derivative2_ξ₀ = 2(d2a_ξ₀ * a_ξ₀ + α_ξ₀^2 + d2b_ξ₀ * b_ξ₀ + β_ξ₀^2)

        !Arblib.contains_zero(abs2_Q_derivative2_ξ₀)
    end

    if verified_zero
        verbose && @info "Verified monotonicity on (0, ξ₀)" ξ₀
    else
        verbose && @warn "Could not verify monotonicity on (0, ξ₀)" ξ₀
        return false, Arb[], Bool[]
    end

    if Arblib.contains_zero(abs2_Q_derivatives[end])
        # We can never verify the existence of a critical point in the
        # last interval since we need to check that the sign is
        # non-zero to the right of it. If the enclosure of the
        # derivative for last interval contains zero we can therefore
        # fail early.
        verbose && @warn "Could not verify monotonicity on last subinterval of [ξ₀, ξ₁]"
        return false, Arb[], Bool[]
    end

    # Count critical points on [ξ₀, ξ₁]
    zeros, verified = let
        # Find all intervals on [ξ₀, ξ₁] for which the enclosure of
        # the derivative contains zero. These are potential critical
        # points.
        zeros = filter(>=(i), findall(Arblib.contains_zero, abs2_Q_derivatives))

        if isempty(zeros)
            Arb[], Bool[]
        else
            # Group the intervals with potential critical points into
            # consecutive chunks

            # End indices for all chunks
            zero_chunks_end_indices = pushfirst!(
                findall(push!((zeros .+ 1)[1:(end-1)] .!= zeros[2:end], true)),
                0,
            )

            # Each chunk as a vector of enclosures
            zero_chunks = map(2:lastindex(zero_chunks_end_indices)) do i
                zeros[zero_chunks_end_indices[i-1]+1]:zeros[zero_chunks_end_indices[i]]
            end

            # Filter out the chunks where the second derivative is
            # non-zero, so the derivative is monotone, and the
            # endpoints have the same sign. There can't be a critical
            # point there.
            zero_chunks = filter(zero_chunks) do zero_chunk
                !(
                    all(!Arblib.contains_zero, abs2_Q_derivative2s[zero_chunk]) &&
                    (
                        Arblib.sgn_nonzero(abs2_Q_derivatives[zero_chunk[1]-1]) ==
                        Arblib.sgn_nonzero(abs2_Q_derivatives[zero_chunk[end]+1])
                    )
                )
            end

            # Union of enclosures for each chunk
            zeros_chunked = map(zero_chunks) do zero_chunk
                reduce(Arblib.union, ξs[zero_chunk])
            end

            verified_zeros = map(zero_chunks) do zero_chunk
                # Check that second derivative is non-zero and that
                # the sign at the endpoints differs
                all(!Arblib.contains_zero, abs2_Q_derivative2s[zero_chunk]) && (
                    Arblib.sgn_nonzero(abs2_Q_derivatives[zero_chunk[1]-1]) *
                    Arblib.sgn_nonzero(abs2_Q_derivatives[zero_chunk[end]+1]) ==
                    -1
                )
            end

            zeros_chunked, verified_zeros
        end
    end

    if all(verified)
        verbose && @info "Verified $(length(zeros)) critical points on [ξ₀, ξ₁]"
    else
        verbose && @warn "Could not verify all critical points on [ξ₀, ξ₁]" zeros verified
    end

    return verified_zero & all(verified), zeros, verified
end
