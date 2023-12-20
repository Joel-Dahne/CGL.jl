G_solve(μ₀, γ₀, κ₀, ξ₁, λ::CGLParams; verbose = false) = G_solve(
    convert(Arb, μ₀),
    convert(Arb, real(γ₀)),
    convert(Arb, imag(γ₀)),
    convert(Arb, κ₀),
    convert(Arb, ξ₁),
    CGLParams{Arb}(λ);
    verbose,
)

function G_solve(
    μ₀::Arb,
    γ₀_real::Arb,
    γ₀_imag::Arb,
    κ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    verbose = false,
)
    x₀ = SVector(μ₀, γ₀_real, γ₀_imag, κ₀)

    # Pick radius as large as possible but so that J \ y can still be
    # computed
    y = ArbMatrix(G_real(x₀..., ξ₁, λ))

    rs = Arb(10) .^ range(-5, -10, 16)

    # Given an r check if J \ y can be computed
    is_ok(r::Arb) =
        let x = add_error.(x₀, r), J = ArbMatrix(G_jacobian_real(x..., ξ₁, λ))
            !iszero(Arblib.solve!(similar(y), J, y))
        end
    # Needed for searchsortedfirst to work correctly. It applies
    # is_ok also to the true argument.
    is_ok(b::Bool) = b

    rs_idx = searchsortedfirst(rs, true, by = is_ok)

    if rs_idx > lastindex(rs)
        verbose && @error "Could not invert with smallest considered radius" rs[end]
        return SVector(
            indeterminate(μ₀),
            indeterminate(μ₀),
            indeterminate(μ₀),
            indeterminate(μ₀),
        )
    end

    if rs_idx == firstindex(rs)
        verbose && @info "Could invert with largest considered radius"
    end

    r = rs[rs_idx]

    verbose && @info "Found optimal radius" r

    # Actual ball we work with
    x = SVector(
        add_error(x₀[1], r),
        add_error(x₀[2], 100r),
        add_error(x₀[3], 100r),
        add_error(x₀[4], r),
    )

    G = x -> CGL.G_real(x..., ξ₁, λ)
    dG = x -> CGL.G_jacobian_real(x..., ξ₁, λ)

    res = CGL.verify_and_refine_root(G, dG, x; verbose)

    return res
end
