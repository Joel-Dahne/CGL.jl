function run_branch_points_verification(directory::AbstractString; verbose = true)
    subdirectories = readdir(directory)
    @assert all(subdir -> isdir(joinpath(directory, subdir)), subdirectories)

    verbose && @info "Found $(length(subdirectories)) subdirectories"

    verbose && @info "Reading subdirectories"

    λs_fix_kappa_all_branches = map(subdirectories) do subdir
        _, λ, fix_kappa = read_parameters(joinpath(directory, subdir, "parameters.csv"))

        verbose && @info "Read from subdirectory $subdir" λ.ω

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

        rescaled_values = map(λs, values) do λ, value
            map(value) do μ_γ_κ
                scale_params(μ_γ_κ..., λ, scaling = inv(sqrt(λ.ω)))
            end
        end

        rescaled_values_overlaps = true
        for k in 1:length(rescaled_values)
            for l in k+1:length(rescaled_values)
                v1 = rescaled_values[k]
                v2 = rescaled_values[l]
                rescaled_values_overlaps |= all(zip(v1, v2)) do (μ_γ_κ_1, μ_γ_κ_2)
                    all(Arblib.overlaps.(μ_γ_κ_1, μ_γ_κ_2))
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
