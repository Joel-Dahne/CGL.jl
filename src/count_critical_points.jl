"""
count_critical_points(μ::Arb, γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}; verbose)

**IMPROVE:** We could let the rightmost interval use a lower bound
which is larger than ξ₁. This would require slightly rewording the
section in the paper. But it avoids having to use a larger ξ₁ when
computing the solution.

**IMPROVE:** The enclosures for the derivatives of abs2(Q) are not
sufficiently good to isolate all critical points for `j >= 5`. One
option for improving the enclosures might be to compute them directly
in CAPD.
"""
function count_critical_points(
    μ::Arb,
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    verbose = false,
)
    # Find ξ₂ such that monotonicity is verified on [ξ₂, ∞)
    ξ₂ = verify_monotonicity_infinity(γ, κ, ϵ, ξ₁, λ; verbose)

    if isfinite(ξ₂)
        verbose && @info "Verified monotonicity on [ξ₂, ∞)" ξ₂
    else
        verbose && @warn "Could not verify monotonicity on [ξ₂, ∞) for any ξ₂"
        return false, Arb[], Bool[]
    end

    if ξ₂ > 2ξ₁
        verbose && @warn "Got ξ₂ > 2ξ₁ - aborting early" ξ₁
        return false, Arb[], Bool[]
    end

    # Compute enclosure on [0, ξ₂]

    ξs, Qs, d2Qs, abs2_Q_derivatives, abs2_Q_derivative2s =
        Q_zero_capd_curve(μ, κ, ϵ, ξ₂, λ)

    if !all(Q -> all(isfinite, Q), Qs)
        verbose && @warn "Could not enclose curve on [0, ξ₂]"
        return false, Arb[], Bool[]
    end

    # Pick ξ₀ value
    i = findfirst(!Arblib.contains_zero, abs2_Q_derivatives)
    ξ₀ = ubound(Arb, ξs[i-1])
    @assert i >= 2

    # Verify monotonicity on [0, ξ₀]
    verified_zero = if λ.d != 3
        # Just check that second order derivative is non-zero
        all(!Arblib.contains_zero, abs2_Q_derivative2s[1:i-1])
    else
        # For d = 3 we get very bad bounds for d2Qs near zero. Instead
        # we evaluate the Taylor expansion directly on [0, ξ₀] to
        # bound it
        (a_ξ₀, b_ξ₀, α_ξ₀, β_ξ₀), (d2a_ξ₀, d2b_ξ₀) =
            Q_zero_taylor(μ, κ, ϵ, ξ₀, λ, enclose_curve = Val{true}())

        abs2_Q_derivative2_ξ₀ = 2(d2a_ξ₀ * a_ξ₀ + α_ξ₀^2 + d2b_ξ₀ * b_ξ₀ + β_ξ₀^2)

        !Arblib.contains_zero(abs2_Q_derivative2_ξ₀)
    end

    if verified_zero
        verbose && @info "Verified monotonicity on [0, ξ₀]" ξ₀
    else
        verbose && @warn "Could not verify monotonicity on [0, ξ₀]" ξ₀
    end

    # Count critical points on [ξ₀, ξ₂]
    zeros, verified = let
        # Find all intervals on [ξ₀, ξ₂] for which the enclosure
        # contains zero
        zeros = filter(>=(i), findall(Arblib.contains_zero, abs2_Q_derivatives))

        if isempty(zeros)
            Arb[], Bool[]
        else
            # Group the intervals into consecutive chunks

            # End indices for all chunks
            zero_chunks_end_indices = pushfirst!(
                findall(push!((zeros.+1)[1:end-1] .!= zeros[2:end], true)),
                0,
            )

            # Each chunk as a vector of enclosures
            zero_chunks = map(2:lastindex(zero_chunks_end_indices)) do i
                zeros[zero_chunks_end_indices[i-1]+1]:zeros[zero_chunks_end_indices[i]]
            end

            # Union of enclosures for each chunk
            zeros_chunked = map(zero_chunks) do zero_chunk
                reduce(Arblib.union, ξs[zero_chunk])
            end

            # IMPROVE: Handle case with false zero and where second
            # derivative can verify not a zero
            verified_zeros = map(zero_chunks) do zero_chunk
                # Check that derivative is non-zero and that the sign
                # at the endpoints differs
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
        verbose && @info "Verified $(length(zeros)) critical points on [ξ₀, ξ₂]"
    else
        verbose && @warn "Could not verify all critical points on [ξ₀, ξ₂]" zeros verified
    end

    return verified_zero & all(verified), zeros, verified
end
