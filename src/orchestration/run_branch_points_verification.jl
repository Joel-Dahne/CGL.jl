function run_branch_points_verification(directory::AbstractString; verbose = true)
    subdirectories = readdir(directory)
    @assert all(subdir -> isdir(joinpath(directory, subdir)), subdirectories)

    verbose && @info "Found $(length(subdirectories)) subdirectories"

    verbose && @info "Reading subdirectories"

    λs_fix_kappa_all_branches = map(subdirectories) do subdir
        _, λ, fix_kappa, ξ₁_strategy, ξ₁_strategy_value =
            read_parameters(joinpath(directory, subdir, "parameters.csv"))

        verbose &&
            @info "Read from subdirectory $subdir" λ.ω ξ₁_strategy ξ₁_strategy_value

        files =
            filter(startswith("branch_points"), readdir(joinpath(directory, subdir)))

        branches = map(files) do file
            read_branch_points_csv(joinpath(directory, subdir, file))
        end

        λ, fix_kappa, branches
    end

    λs = getindex.(λs_fix_kappa_all_branches, 1)
    fix_kappas = getindex.(λs_fix_kappa_all_branches, 2)
    all_branches = getindex.(λs_fix_kappa_all_branches, 3)

    verbose && @info "Checking that all subdirectories contain compatible data"

    if !allequal(λ -> (λ.d, λ.δ, λ.σ), λs)
        # IMPROVE: Should we allow recovering from this state?
        error("not all subdirectories have the same values for d, δ and σ")
    end

    if !allequal(length.(all_branches))
        # IMPROVE: Should we allow recovering from this state?
        error("not all subdirectories have the same number of branches")
    end

    if !allequal(map(branches -> nrow.(branches), all_branches))
        # IMPROVE: Should we allow recovering from this state?
        error("not all subdirectories have the same length for the branches")
    end

    # TODO: Split into branches which fixes ϵ and branches which fixes κ
    @assert all(!, fix_kappas) # FIXME: Handle other cases

    verbose && @info "Checking the results for all subdirectories overlap"

    # TODO: Check that enclosures of branches agree
    for j in eachindex(all_branches[1])
        branches_j = getindex.(all_branches, j)

        # Check that they all have the same ϵ values
        # TODO: Handle other cases?
        @assert allequal(branch -> branch.ϵ, branches_j)

        values = map(branch -> tuple.(branch.μ, branch.γ, branch.κ), branches_j)
        ξ₁s = map(branch -> branch.ξ₁, branches_j)

        rescaled_values = map(λs, values) do λ, value
            map(value) do μ_γ_κ
                scale_params(μ_γ_κ..., λ, scaling = inv(sqrt(λ.ω)))
            end
        end
        rescaled_ξ₁s = map((ξ₁, λ) -> ξ₁ * sqrt(λ.ω), ξ₁s, λs)

        rescaled_values_overlaps = true
        for k = 1:length(branches_j)
            for l = k+1:length(branches_j)
                for (v_1, v_2, ξ₁_1, ξ₁_2) in zip(
                    rescaled_values[k],
                    rescaled_values[l],
                    rescaled_ξ₁s[k],
                    rescaled_ξ₁s[l],
                )

                    μ1, γ1, κ1 = v_1
                    μ2, γ2, κ2 = v_2

                    if Arblib.overlaps(ξ₁_1, ξ₁_2)
                        overlaps = all(Arblib.overlaps.(v_1, v_2))
                    else
                        overlaps = Arblib.overlaps(μ1, μ2) && Arblib.overlaps(κ1, κ2)
                    end

                    if !overlaps
                        @error "Subdirectories $k and $l don't agree" v_1 v_2 Arblib.overlaps.(
                            v_1,
                            v_2,
                        )
                        rescaled_values_overlaps = false
                    end
                end
            end
        end

        if rescaled_values_overlaps
            verbose && @info "All results overlap for j = $j"
        else
            @error "Results don't overlap for j = $j"
        end
    end

    return all_branches
end
