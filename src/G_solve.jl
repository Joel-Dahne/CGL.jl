G_solve(μ₀, γ₀, κ₀, ϵ, ξ₁, λ::CGLParams; rs = 10 .^ range(-5, -10, 8), verbose = false) =
    G_solve(
        convert(Arb, μ₀),
        convert(Arb, real(γ₀)),
        convert(Arb, imag(γ₀)),
        convert(Arb, κ₀),
        convert(Arb, ϵ),
        convert(Arb, ξ₁),
        CGLParams{Arb}(λ),
        rs = convert(Vector{Arb}, rs);
        verbose,
    )

# Work in progress
function G_solve_alt(
    μ₀::Arb,
    γ₀_real::Arb,
    γ₀_imag::Arb,
    κ₀::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    max_iterations::Integer = 10,
    verbose = false,
    return_uniqueness::Union{Val{false},Val{true}} = Val{false}(),
)
    @assert return_uniqueness isa Val{false}

    x = SVector(μ₀, γ₀_real, γ₀_imag, κ₀)

    if any(!isfinite, x)
        verbose && @error "Non-finite input"
        return indeterminate.(x)
    end

    verbose && @info "Iteration 0" x

    G_x = x -> G(x..., ϵ, ξ₁, λ)
    dG_x = x -> G_jacobian_kappa(x..., ϵ, ξ₁, λ)

    isproved = false
    for i = 1:max_iterations
        x_new = newton_step(G_x, dG_x, x)

        if !all(isfinite, x_new)
            verbose && @warn "Newton step failed" x_new
            break
        end

        if all(Arblib.contains_interior.(x, x_new))
            verbose && @info "Success" x_new
            x = x_new
            isproved = true
            break
        end

        verbose && @info "Iteration $i" x_new

        x = add_error.(x_new, 0.01radius.(x_new))
    end

    if isproved
        return x
    else
        verbose && @warn "Could not prove root"
        return indeterminate.(x)
    end
end


# Solve with ϵ fixed
function G_solve(
    μ₀::Arb,
    γ₀_real::Arb,
    γ₀_imag::Arb,
    κ₀::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    rs::Vector{Arb} = Arb(10) .^ range(-5, -10, 8),
    try_double_radius = true,
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
    y₀ = ArbMatrix(G(x₀..., ϵ, ξ₁, λ))

    # Scale the arguments according to the norms of the columns of the
    # Jacobian. Taking such that the scale for κ is 1
    J₀ = G_jacobian_kappa(x₀..., ϵ, ξ₁, λ)

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

# Solve with κ fixed
function G_solve_fix_kappa(
    μ₀::Arb,
    γ₀_real::Arb,
    γ₀_imag::Arb,
    κ::Arb,
    ϵ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    rs::Vector{Arb} = Arb(10) .^ range(-5, -10, 8),
    try_double_radius = true,
    verbose = false,
    return_uniqueness::Union{Val{false},Val{true}} = Val{false}(),
)
    x₀ = SVector(μ₀, γ₀_real, γ₀_imag, ϵ₀)

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
    y₀ = ArbMatrix(G(x₀[1:3]..., κ, x₀[4], ξ₁, λ))

    # Scale the arguments according to the norms of the columns of the
    # Jacobian. Taking such that the scale for ϵ is 1
    J₀ = G_jacobian_epsilon(x₀[1:3]..., κ, x₀[4], ξ₁, λ)

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
