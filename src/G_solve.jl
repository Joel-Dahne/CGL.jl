function G_solve_fix_epsilon_alt(
    μ::Arb,
    γ_real::Arb,
    γ_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    return_uniqueness::Union{Val{false},Val{true}} = Val{false}(),
    try_expand_uniqueness = return_uniqueness isa Val{true},
    expansion_rate = 0.05,
    max_iterations = 10,
    verbose = false,
    extra_verbose = false,
)
    G_x = x -> G(x..., ϵ, ξ₁, λ)
    dG_x = x -> G_jacobian_kappa(x..., ϵ, ξ₁, λ)

    root, root_uniqueness = verify_root_from_approximation(
        G_x,
        dG_x,
        SVector(μ, γ_real, γ_imag, κ);
        expansion_rate,
        max_iterations,
        verbose,
    )

    if try_expand_uniqueness
        verbose && @info "Expanding region for uniqueness"
        root_uniqueness =
            expand_uniqueness(G_x, dG_x, root_uniqueness; verbose, extra_verbose)
    end

    if return_uniqueness isa Val{true}
        return root, root_uniqueness
    else
        return root
    end
end

function G_solve_fix_kappa_alt(
    μ::Arb,
    γ_real::Arb,
    γ_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    return_uniqueness::Union{Val{false},Val{true}} = Val{false}(),
    try_expand_uniqueness = return_uniqueness isa Val{true},
    expansion_rate = 0.05,
    max_iterations::Integer = 10,
    verbose = false,
    extra_verbose = false,
)
    G_x = x -> G(x[1:3]..., κ, x[4], ξ₁, λ)
    dG_x = x -> G_jacobian_epsilon(x[1:3]..., κ, x[4], ξ₁, λ)

    root, root_uniqueness = verify_root_from_approximation(
        G_x,
        dG_x,
        SVector(μ, γ_real, γ_imag, ϵ);
        expansion_rate,
        max_iterations,
        verbose,
    )

    if try_expand_uniqueness
        verbose && @info "Expanding region for uniqueness"
        root_uniqueness =
            expand_uniqueness(G_x, dG_x, root_uniqueness; verbose, extra_verbose)
    end

    if return_uniqueness isa Val{true}
        return root, root_uniqueness
    else
        return root
    end
end
