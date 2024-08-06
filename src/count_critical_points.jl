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
    ξs, Qs, d2Qs = Q_zero_capd_curve(μ, κ, ϵ, ξ₁, λ)

    # Compute enclosures of the first and second derivative of abs2(Q)
    # (divided by 2)
    abs2_Q_derivatives = map(Qs) do (a, b, α, β)
        α * a + β * b
    end

    abs2_Q_derivative2s = map(ξs, Qs, d2Qs) do ξ, (a, b, α, β), (d2a, d2b)
        if λ.d == 3
            # For d = 3 the enclosures for d2a and d2b are sometimes quite
            # bad, unclear why. To get better enclosures we intersect with
            # the value we get by directly evaluating the ODE
            _, _, d2a_alt, d2b_alt = cgl_equation_real(SVector(a, b, α, β), κ, ϵ, ξ, λ)
            @assert Arblib.overlaps(d2a, d2a_alt)
            @assert Arblib.overlaps(d2b, d2b_alt)
            if isfinite(d2a_alt)
                d2a = Arblib.intersection(d2a, d2a_alt)
            end
            if isfinite(d2b_alt)
                d2b = Arblib.intersection(d2b, d2b_alt)
            end
        end

        d2a * a + α^2 + d2b * b + β^2
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

        abs2_Q_derivative2_ξ₀ = d2a_ξ₀ * a_ξ₀ + α_ξ₀^2 + d2b_ξ₀ * b_ξ₀ + β_ξ₀^2

        !Arblib.contains_zero(abs2_Q_derivative2_ξ₀)
    end

    if verified_zero
        verbose && @info "Verified monotonicity on [0, ξ₀]" ξ₀
    else
        verbose && @warn "Could not verify monotonicity on [0, ξ₀]" ξ₀
    end

    # Verify monotonicity on [ξ₁, ∞)
    verified_infinity = verify_monotonicity_infinity(γ, κ, ϵ, ξ₁, λ; verbose)

    if verified_infinity
        verbose && @info "Verified monotonicity on [ξ₁, ∞)"
    else
        verbose && @warn "Could not verify monotonicity on [ξ₁, ∞)"
    end

    # Count critical points on [ξ₀, ξ₁]
    zeros, verified = let
        # Find all intervals on [ξ₀, ξ₁] for which the enclosure
        # contains zero
        zeros = filter(>=(i), findall(Arblib.contains_zero, abs2_Q_derivatives))

        if isempty(zeros)
            Arb[], Bool[]
        else
            # Group the intervals into consecutive chunks
            zero_chunks_end_indices = pushfirst!(
                findall(push!((zeros.+1)[1:end-1] .!= zeros[2:end], true)),
                0,
            )

            zero_chunks = map(2:lastindex(zero_chunks_end_indices)) do i
                zeros[zero_chunks_end_indices[i-1]+1]:zeros[zero_chunks_end_indices[i]]
            end

            zeros_chunked = map(zero_chunks) do zero_chunk
                reduce(Arblib.union, ξs[zero_chunk])
            end

            # IMPROVE: Handle case with false zero and where second
            # derivative can verify not a zero
            verified_zeros = map(zero_chunks) do zero_chunk
                # Check that derivative is non-zero and that the sign
                # at the endpoints differs
                all(!Arblib.contains_zero, abs2_Q_derivative2s[zero_chunk]) &&
                    Arblib.sgn_nonzero(abs2_Q_derivatives[zero_chunk[1]-1]) *
                    Arblib.sgn_nonzero(abs2_Q_derivatives[zero_chunk[end]+1]) == -1
            end

            zeros_chunked, verified_zeros
        end
    end

    if all(verified)
        verbose && @info "Verified $(length(zeros)) critical points on [ξ₀, ξ₁]"
    else
        verbose && @warn "Could not verify all critical points on [ξ₀, ξ₁]" zeros verified
    end

    return verified_zero & verified_infinity & all(verified), zeros, verified
end
