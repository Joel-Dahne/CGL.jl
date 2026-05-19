### A Pluto.jl notebook ###
# v0.20.25

using Markdown
using InteractiveUtils

# ╔═╡ 031b0a42-9f10-11ee-28d2-5957efedccab
begin
    using Pkg
    Pkg.activate("..", io = devnull)
    using Arblib
    using CGL
    using DataFrames
    using LaTeXStrings
    using OhMyThreads
    using Plots
    using PlutoUI

    setprecision(Arb, 128)

    nothing
end

# ╔═╡ c9556190-f0c5-44f8-a967-e27bd80829db
md"""
# Results for the NLS equation

This notebook generates the results for the NLS equation. It computes all the required data and generates the tables given in Section 3.
"""

# ╔═╡ b7903f77-8cd1-455f-aec7-be142ce71c4e
TableOfContents()

# ╔═╡ 6ea53a4c-7365-4599-b4e4-a4633c31063d
md"""
## Case I

For Case I we prove the existence of 8 self-similar profiles.
"""

# ╔═╡ afe77be4-3366-45b5-b3eb-a1413c12c3d6
parameters_1 = [
    (j = 1, d = 1),
    (j = 2, d = 1),
    (j = 3, d = 1),
    (j = 4, d = 1),
    (j = 5, d = 1),
    (j = 6, d = 1),
    (j = 7, d = 1),
    (j = 8, d = 1),
]

# ╔═╡ 6ef41ef4-cf41-4e2e-83f3-8f86afb6c120
md"""
Prove the existence of a self-similar profile for each set of parameters.
"""

# ╔═╡ a0609434-aeea-40f8-85e0-36ff26046aea
res_1 = OhMyThreads.tmap(parameters_1) do (j, d)
    μ, γ, κ, ϵ, ξ₁, Λ = CGL.sverak_params(Arb, j, d)
    CGL.G_solve_fix_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁, Λ), ξ₁
end

# ╔═╡ 2d1d2b8e-775f-4722-aa73-7f6e8a18ded5
md"""
Verify that all of the computations were successful.
"""

# ╔═╡ acccef27-6b2e-4e57-8020-4d0c3815aa79
@assert_proof all(root -> all(isfinite, root), getindex.(res_1, 1))

# ╔═╡ fb8bd43a-efde-4b8d-9968-0d315e631331
md"""
Count the number of critical points for each self-similar profile.
"""

# ╔═╡ f7d24e88-fe9a-4ff7-98dd-6db0d9d95736
critical_points_1 = OhMyThreads.tmap(res_1, parameters_1) do (sol, ξ₁), (j, d)
    μ, γ_real, γ_imag, κ = sol
    ϵ = zero(μ)
    Λ = CGL.sverak_params(Arb, j, d)[6]
    CGL.count_critical_points(μ, Acb(γ_real, γ_imag), κ, ϵ, ξ₁, Λ)
end

# ╔═╡ 5b01f78d-1986-494e-9bc8-77a8e6f4eeac
md"""
Verify that the number of critical points were computed for all profiles.
"""

# ╔═╡ 820cf8ac-cc01-48c5-8bfa-5e8ae9d0954d
@assert_proof all(getindex.(critical_points_1, 1))

# ╔═╡ 169728ed-a34c-4350-a4db-23e04c331deb
md"""
In the paper we say that the number of critical points is given by ``j - 1``. Verify that this is the case.
"""

# ╔═╡ 0143b3b3-a1fb-463e-bcf5-70790b651ea0
@assert_proof length.(getindex.(critical_points_1, 2)) == getindex.(parameters_1, 1) .- 1

# ╔═╡ 19a76800-f3a6-4541-89d9-90d31c04b32e
md"""
We collect all data in a dataframe to make it easier to read.
"""

# ╔═╡ e2b0a9b6-31b3-443e-b1f2-872d0739b904
df_1 = let df = DataFrame()
    df.j = getfield.(parameters_1, :j)
    df.d = getfield.(parameters_1, :d)
    df.μ = getindex.(getindex.(res_1, 1), 1)
    df.γ_real = getindex.(getindex.(res_1, 1), 2)
    df.γ_imag = getindex.(getindex.(res_1, 1), 3)
    df.κ = getindex.(getindex.(res_1, 1), 4)
    df.ξ₁ = Float64.(getindex.(res_1, 2))
    df.num_critical_points = ifelse.(
        getindex.(critical_points_1, 1),
        length.(getindex.(critical_points_1, 2)),
        missing,
    )
    df
end

# ╔═╡ 285ba40d-98bd-4c23-9cd6-74fd3c5591df
md"""
## Case II

For Case II we prove the existence of 2 self-similar profiles.
"""

# ╔═╡ 54c46505-cf47-404c-aa57-f331724951ee
parameters_2 = [(j = 1, d = 3), (j = 2, d = 3)]

# ╔═╡ 38afdbd4-777f-45c2-9a1c-add9b755caf5
md"""
Prove the existence of a self-similar profile for each set of parameters.
"""

# ╔═╡ d836fdb5-8ea8-44a7-9c57-f4656545b27b
res_2 = OhMyThreads.tmap(parameters_2) do (j, d)
    μ, γ, κ, ϵ, ξ₁, Λ = CGL.sverak_params(Arb, j, d)
    CGL.G_solve_fix_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁, Λ), ξ₁
end

# ╔═╡ c51e3358-22f6-4ba6-b720-3d79f6311e7d
md"""
Verify that all of the computations were successful.
"""

# ╔═╡ 81de63ac-8239-4b74-9b1b-4573323c0a8d
@assert_proof all(root -> all(isfinite, root), getindex.(res_2, 1))

# ╔═╡ 6eff2b62-7b70-404c-a25c-2ce660b8d454
md"""
Count the number of critical points for each self-similar profile.
"""

# ╔═╡ cab3f05d-c949-490c-9367-6926732df899
critical_points_2 = OhMyThreads.tmap(res_2, parameters_2) do (sol, ξ₁), (j, d)
    μ, γ_real, γ_imag, κ = sol
    ϵ = zero(μ)
    Λ = CGL.sverak_params(Arb, j, d)[6]
    CGL.count_critical_points(μ, Acb(γ_real, γ_imag), κ, ϵ, ξ₁, Λ)
end

# ╔═╡ ce85bc5b-cc2d-44fd-82b5-67c4cf568a6d
md"""
Verify that the number of critical points were computed for all profiles.
"""

# ╔═╡ 170116be-d090-4189-a700-49ae2e5c5cc8
@assert_proof all(getindex.(critical_points_2, 1))

# ╔═╡ 519fcf2e-cf20-400a-8c61-3bd425de01d8
md"""
In the paper we say that the number of critical points is given by ``j - 1``. Verify that this is the case.
"""

# ╔═╡ f3f89f57-0ebe-4685-992c-78da9e47dba2
@assert_proof length.(getindex.(critical_points_2, 2)) == getindex.(parameters_2, 1) .- 1

# ╔═╡ 214200b4-6cd2-4177-b61a-40ace25cd58f
md"""
We collect all data in a dataframe to make it easier to read.
"""

# ╔═╡ 2e75ba82-dbe8-472d-b7ab-e77175b1699a
df_2 = let df = DataFrame()
    df.j = getfield.(parameters_2, :j)
    df.d = getfield.(parameters_2, :d)
    df.μ = getindex.(getindex.(res_2, 1), 1)
    df.γ_real = getindex.(getindex.(res_2, 1), 2)
    df.γ_imag = getindex.(getindex.(res_2, 1), 3)
    df.κ = getindex.(getindex.(res_2, 1), 4)
    df.ξ₁ = Float64.(getindex.(res_2, 2))
    df.num_critical_points = ifelse.(
        getindex.(critical_points_2, 1),
        length.(getindex.(critical_points_2, 2)),
        missing,
    )
    df
end

# ╔═╡ 947a1fc3-6d56-4349-9cfa-1096a32f9241
md"""
## LaTeX output

Prepare formatted LaTeX output for the paper. These are the two tables in Section 3.
"""

# ╔═╡ eecc2ec2-ffb1-434a-a6a8-490672b56284
df_1_string = let df = DataFrame()
    df.j = df_1.j
    df.μ = CGL.format_interval_precise.(df_1.μ)
    df.γ = CGL.format_interval_precise.(Acb.(df_1.γ_real, df_1.γ_imag))
    df.κ = CGL.format_interval_precise.(df_1.κ)
    df.ξ₁ = Int.(df_1.ξ₁)
    df
end

# ╔═╡ a2449e4c-87f9-44ef-a1af-0eac9581f027
df_2_string = let df = DataFrame()
    df.j = df_2.j
    df.μ = CGL.format_interval_precise.(df_2.μ)
    df.γ = CGL.format_interval_precise.(Acb.(df_2.γ_real, df_2.γ_imag))
    df.κ = CGL.format_interval_precise.(df_2.κ)
    df.ξ₁ = Int.(df_2.ξ₁)
    df
end

# ╔═╡ 36e2d4b8-60bd-418c-8311-328b1f6ed937
map(eachrow(df_1_string)) do (j, μ, γ, κ, ξ₁)
    "\\($j\\) & \\($μ\\) & \\($γ\\) & \\($κ\\) & \\($ξ₁\\)\\\\\n"
end |> join |> println

# ╔═╡ e9d7b985-3e7c-4900-9d36-511d07ee6d6c
map(eachrow(df_2_string)) do (j, μ, γ, κ, ξ₁)
    "\\($j\\) & \\($μ\\) & \\($γ\\) & \\($κ\\) & \\($ξ₁\\)\\\\\n"
end |> join |> println

# ╔═╡ Cell order:
# ╟─c9556190-f0c5-44f8-a967-e27bd80829db
# ╠═031b0a42-9f10-11ee-28d2-5957efedccab
# ╠═b7903f77-8cd1-455f-aec7-be142ce71c4e
# ╟─6ea53a4c-7365-4599-b4e4-a4633c31063d
# ╠═afe77be4-3366-45b5-b3eb-a1413c12c3d6
# ╟─6ef41ef4-cf41-4e2e-83f3-8f86afb6c120
# ╠═a0609434-aeea-40f8-85e0-36ff26046aea
# ╟─2d1d2b8e-775f-4722-aa73-7f6e8a18ded5
# ╠═acccef27-6b2e-4e57-8020-4d0c3815aa79
# ╟─fb8bd43a-efde-4b8d-9968-0d315e631331
# ╠═f7d24e88-fe9a-4ff7-98dd-6db0d9d95736
# ╟─5b01f78d-1986-494e-9bc8-77a8e6f4eeac
# ╠═820cf8ac-cc01-48c5-8bfa-5e8ae9d0954d
# ╟─169728ed-a34c-4350-a4db-23e04c331deb
# ╠═0143b3b3-a1fb-463e-bcf5-70790b651ea0
# ╟─19a76800-f3a6-4541-89d9-90d31c04b32e
# ╠═e2b0a9b6-31b3-443e-b1f2-872d0739b904
# ╟─285ba40d-98bd-4c23-9cd6-74fd3c5591df
# ╠═54c46505-cf47-404c-aa57-f331724951ee
# ╟─38afdbd4-777f-45c2-9a1c-add9b755caf5
# ╠═d836fdb5-8ea8-44a7-9c57-f4656545b27b
# ╟─c51e3358-22f6-4ba6-b720-3d79f6311e7d
# ╠═81de63ac-8239-4b74-9b1b-4573323c0a8d
# ╟─6eff2b62-7b70-404c-a25c-2ce660b8d454
# ╠═cab3f05d-c949-490c-9367-6926732df899
# ╟─ce85bc5b-cc2d-44fd-82b5-67c4cf568a6d
# ╠═170116be-d090-4189-a700-49ae2e5c5cc8
# ╟─519fcf2e-cf20-400a-8c61-3bd425de01d8
# ╠═f3f89f57-0ebe-4685-992c-78da9e47dba2
# ╟─214200b4-6cd2-4177-b61a-40ace25cd58f
# ╠═2e75ba82-dbe8-472d-b7ab-e77175b1699a
# ╟─947a1fc3-6d56-4349-9cfa-1096a32f9241
# ╟─eecc2ec2-ffb1-434a-a6a8-490672b56284
# ╟─a2449e4c-87f9-44ef-a1af-0eac9581f027
# ╠═36e2d4b8-60bd-418c-8311-328b1f6ed937
# ╠═e9d7b985-3e7c-4900-9d36-511d07ee6d6c
