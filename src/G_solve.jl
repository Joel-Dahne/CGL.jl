function G_solve_fix_epsilon(
    μ::Arb,
    γ_real::Arb,
    γ_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    rs::Vector{Arb} = Arb(10) .^ range(-5, -10, 8),
    try_double_radius = true,
    verbose = false,
    return_uniqueness::Union{Val{false},Val{true}} = Val{false}(),
)
    x₀ = SVector(μ, γ_real, γ_imag, κ)

    if any(!isfinite, x₀)
        verbose && @error "Non-finite input"
        res = SVector(
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
        )
        if return_uniqueness isa Val{false}
            return res
        else
            return res, deepcopy(res)
        end
    end

    # Pick radius as large as possible but so that J \ y can still be
    # computed
    y₀ = ArbMatrix(G(x₀..., ϵ, ξ₁, λ))

    # Scale the arguments according to the norms of the columns of the
    # Jacobian. Taking such that the scale for κ is 1
    J₀ = G_jacobian_kappa(x₀..., ϵ, ξ₁, λ)

    if iszero(Arblib.solve!(similar(y₀), ArbMatrix(J₀), y₀))
        verbose && @error "Could not invert with zero radius"
        res = SVector(
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
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
        let x = add_error.(x₀, r * r_scaling),
            J = ArbMatrix((G_jacobian_kappa(x..., ϵ, ξ₁, λ)))

            !iszero(Arblib.solve!(similar(y₀), J, y₀))
        end
    # Needed for searchsortedfirst to work correctly. It applies
    # is_ok also to the true argument.
    is_ok(b::Bool) = b

    rs_idx = searchsortedfirst(rs, true, by = is_ok)

    if rs_idx > lastindex(rs)
        verbose && @error "Could not invert with smallest considered radius" rs[end]
        res = SVector(
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
        )
        if return_uniqueness isa Val{false}
            return res
        else
            return res, deepcopy(res)
        end
    end

    r = rs[rs_idx]

    if rs_idx == firstindex(rs)
        verbose && @info "Could invert with largest considered radius"
        if try_double_radius
            r_new = 2r
            doubles = 0
            while is_ok(r_new)
                r = r_new
                r_new = 2r
                doubles += 1
            end
            verbose && @info "Doubled radius $doubles times"
        end
    end

    verbose && @info "Found optimal radius" r

    # Actual ball we work with
    x = add_error.(x₀, r * r_scaling)

    G_x = x -> G(x..., ϵ, ξ₁, λ)
    dG_x = x -> G_jacobian_kappa(x..., ϵ, ξ₁, λ)

    res = verify_and_refine_root(G_x, dG_x, x; verbose)

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

# Solve with κ fixed
function G_solve_fix_kappa(
    μ::Arb,
    γ_real::Arb,
    γ_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    rs::Vector{Arb} = Arb(10) .^ range(-5, -10, 8),
    try_double_radius = true,
    verbose = false,
    return_uniqueness::Union{Val{false},Val{true}} = Val{false}(),
)
    x₀ = SVector(μ, γ_real, γ_imag, ϵ)

    if any(!isfinite, x₀)
        verbose && @error "Non-finite input"
        res = SVector(
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
        )
        if return_uniqueness isa Val{false}
            return res
        else
            return res, deepcopy(res)
        end
    end

    # Pick radius as large as possible but so that J \ y can still be
    # computed
    y₀ = ArbMatrix(G(x₀[1:3]..., κ, x₀[4], ξ₁, λ))

    # Scale the arguments according to the norms of the columns of the
    # Jacobian. Taking such that the scale for ϵ is 1
    J₀ = G_jacobian_epsilon(x₀[1:3]..., κ, x₀[4], ξ₁, λ)

    if iszero(Arblib.solve!(similar(y₀), ArbMatrix(J₀), y₀))
        verbose && @error "Could not invert with zero radius"
        res = SVector(
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
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
        let x = add_error.(x₀, r * r_scaling),
            J = ArbMatrix((G_jacobian_epsilon(x[1:3]..., κ, x[4], ξ₁, λ)))

            !iszero(Arblib.solve!(similar(y₀), J, y₀))
        end
    # Needed for searchsortedfirst to work correctly. It applies
    # is_ok also to the true argument.
    is_ok(b::Bool) = b

    rs_idx = searchsortedfirst(rs, true, by = is_ok)

    if rs_idx > lastindex(rs)
        verbose && @error "Could not invert with smallest considered radius" rs[end]
        res = SVector(
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
            indeterminate(μ),
        )
        if return_uniqueness isa Val{false}
            return res
        else
            return res, deepcopy(res)
        end
    end

    r = rs[rs_idx]

    if rs_idx == firstindex(rs)
        verbose && @info "Could invert with largest considered radius"
        if try_double_radius
            r_new = 2r
            doubles = 0
            while is_ok(r_new)
                r = r_new
                r_new = 2r
                doubles += 1
            end
            verbose && @info "Doubled $doubles times"
        end
    end

    verbose && @info "Found optimal radius" r

    # Actual ball we work with
    x = add_error.(x₀, r * r_scaling)

    G_x = x -> G(x[1:3]..., κ, x[4], ξ₁, λ)
    dG_x = x -> G_jacobian_epsilon(x[1:3]..., κ, x[4], ξ₁, λ)

    res = verify_and_refine_root(G_x, dG_x, x; verbose)

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
