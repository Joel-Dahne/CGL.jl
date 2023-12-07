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

    # Pick radius as large as possible but so that the Jacobian is
    # still invertible
    y = ArbMatrix(G_real(x₀..., ξ₁, λ))

    r = Arb(10)^-10

    # Check if this radius works
    x = add_error.(x₀, r)
    J = ArbMatrix(G_jacobian_real(x..., ξ₁, λ))
    success = !iszero(Arblib.solve!(similar(y), J, y))
    if !success
        verbose && @error "Could not invert with initial radius" r J
        return SVector(
            indeterminate(μ₀),
            indeterminate(μ₀),
            indeterminate(μ₀),
            indeterminate(μ₀),
        )
    end

    # Try larger radius
    factor = 16
    for iter = 1:50 # Perform at most 100 steps
        x = add_error.(x₀, r)
        J = ArbMatrix(G_jacobian_real(x..., ξ₁, λ))
        success = !iszero(Arblib.solve!(similar(y), J, y))

        if !success
            r = r / factor
            factor = factor ÷ 2
            factor == 1 && break
        end

        r = factor * r
    end

    verbose && @info "Found optimal radius" r

    # Actual ball we work with
    x = add_error.(x₀, r)

    G = x -> CGL.G_real(x..., ξ₁, λ)
    dG = x -> CGL.G_jacobian_real(x..., ξ₁, λ)

    res = CGL.verify_and_refine_root(G, dG, x; verbose)

    return res
end
