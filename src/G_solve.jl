G_solve(μ₀, γ₀, κ₀, ξ₁, λ::CGLParams; rs = 10 .^ range(-5, -10, 8), verbose = false) =
    G_solve(
        convert(Arb, μ₀),
        convert(Arb, real(γ₀)),
        convert(Arb, imag(γ₀)),
        convert(Arb, κ₀),
        convert(Arb, ξ₁),
        CGLParams{Arb}(λ),
        rs = convert(Vector{Arb}, rs);
        verbose,
    )

function G_solve(
    μ₀::Arb,
    γ₀_real::Arb,
    γ₀_imag::Arb,
    κ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    rs::Vector{Arb} = Arb(10) .^ range(-5, -10, 8),
    verbose = false,
    return_uniqueness::Union{Val{false},Val{true}} = Val{false}(),
)
    x₀ = SVector(μ₀, γ₀_real, γ₀_imag, κ₀)

    if any(!isfinite, x₀)
        verbose && @error "Non-finite input"
        res = SVector(
            indeterminate(μ₀),
            indeterminate(μ₀),
            indeterminate(μ₀),
            indeterminate(μ₀),
        )
        if return_uniqueness isa Val{false}
            return res
        else
            return res, deepcopy(res)
        end
    end

    # Pick radius as large as possible but so that J \ y can still be
    # computed
    y₀ = ArbMatrix(G_real(x₀..., ξ₁, λ))

    # Scale the arguments according to the norms of the columns of the
    # Jacobian. Taking such that the scale for κ is 1
    J₀ = G_jacobian_real(x₀..., ξ₁, λ)

    if iszero(Arblib.solve!(similar(y₀), ArbMatrix(J₀), y₀))
        verbose && @error "Could not invert with zero radius"
        res = SVector(
            indeterminate(μ₀),
            indeterminate(μ₀),
            indeterminate(μ₀),
            indeterminate(μ₀),
        )
        if return_uniqueness isa Val{false}
            return res
        else
            return res, deepcopy(res)
        end
    end

    n1, n2, n3, n4 = norm.(eachcol(J₀))
    r_scaling = n4 ./ SVector(n1, n2, n3, n4)

    # Given an r check if J \ y can be computed
    is_ok(r::Arb) =
        let x = add_error.(x₀, r * r_scaling), J = ArbMatrix((G_jacobian_real(x..., ξ₁, λ)))
            !iszero(Arblib.solve!(similar(y₀), J, y₀))
        end
    # Needed for searchsortedfirst to work correctly. It applies
    # is_ok also to the true argument.
    is_ok(b::Bool) = b

    rs_idx = searchsortedfirst(rs, true, by = is_ok)

    if rs_idx > lastindex(rs)
        verbose && @error "Could not invert with smallest considered radius" rs[end]
        res = SVector(
            indeterminate(μ₀),
            indeterminate(μ₀),
            indeterminate(μ₀),
            indeterminate(μ₀),
        )
        if return_uniqueness isa Val{false}
            return res
        else
            return res, deepcopy(res)
        end
    end

    if rs_idx == firstindex(rs)
        verbose && @info "Could invert with largest considered radius"
    end

    r = rs[rs_idx]

    verbose && @info "Found optimal radius" r

    # Actual ball we work with
    x = add_error.(x₀, r * r_scaling)

    G = x -> G_real(x..., ξ₁, λ)
    dG = x -> G_jacobian_real(x..., ξ₁, λ)

    res = verify_and_refine_root(G, dG, x; verbose)

    if return_uniqueness isa Val{false}
        return res
    elseif all(isfinite, res)
        return res, x
    else
        indeterminate_vector = SVector(
            indeterminate(Arb),
            indeterminate(Arb),
            indeterminate(Arb),
            indeterminate(Arb),
        )
        return res, indeterminate_vector
    end
end
