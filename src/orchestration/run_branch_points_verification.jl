function run_branch_points_verification(
    versions::Vector{Vector{DataFrame}},
    λs::Vector{CGLParams{Arb}};
    verbose = true,
    extra_verbose = false,
)
    verbose && @info "Checking that all data is compatible"

    if !allequal(λ -> (λ.d, λ.δ, λ.σ), λs)
        error("not all λs have the same values for d, δ and σ")
    end

    if !allequal(length.(versions))
        error("not all versions have the same number of branches")
    end

    if !allequal(map(version -> nrow.(version), versions))
        error("not all versions have the same lengths for branches")
    end

    verbose && @info "Checking that results for all versions overlap"

    enclosures_overlaps = map(eachindex(versions[1])) do j
        branches_j = getindex.(versions, j)

        if !allequal(branch -> branch.ϵ, branches_j)
            error("ϵ values don't agree for branches j = $j")
        end

        enclosures = map(branch -> tuple.(branch.μ, branch.γ, branch.κ), branches_j)
        ξ₁s = map(branch -> branch.ξ₁, branches_j)

        rescaled_enclosures = map(λs, enclosures) do λ, enclosure
            map(enclosure) do μ_γ_κ
                scale_params(μ_γ_κ..., λ, scaling = inv(sqrt(λ.ω)))
            end
        end
        rescaled_ξ₁s = map((ξ₁, λ) -> ξ₁ * sqrt(λ.ω), ξ₁s, λs)

        enclosures_overlaps_j = true
        for k = 1:length(branches_j)
            for l = k+1:length(branches_j)
                for (enclosure_k, enclosure_l, ξ₁_k, ξ₁_l) in zip(
                    rescaled_enclosures[k],
                    rescaled_enclosures[l],
                    rescaled_ξ₁s[k],
                    rescaled_ξ₁s[l],
                )

                    if Arblib.overlaps(ξ₁_k, ξ₁_l)
                        overlaps = all(Arblib.overlaps.(enclosure_k, enclosure_l))
                    else
                        μ_k, γ_k, κ_k = enclosure_k
                        μ_l, γ_l, κ_l = enclosure_l
                        overlaps =
                            Arblib.overlaps(μ_k, μ_l) && Arblib.overlaps(κ_k, κ_l)
                    end

                    if !overlaps
                        @error "Versions $k and $l don't agree" enclosure_k enclosure_l Arblib.overlaps.(
                            enclosure_k,
                            enclosure_l,
                        )
                        enclosures_overlaps_j = false
                    end
                end
            end
        end

        num_points = nrow(branches_j[1])
        num_finite_points_per_version = map(branches_j) do branch
            count(isfinite, branch.μ)
        end
        verbose && @info "Checking branch j = $j" num_points

        extra_verbose && @info "Finite points per version" num_finite_points_per_version

        if !enclosures_overlaps_j
            @error "Results don't overlap for j = $j"
        end

        enclosures_overlaps_j
    end

    return enclosures_overlaps
end
