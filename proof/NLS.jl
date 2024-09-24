### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 031b0a42-9f10-11ee-28d2-5957efedccab
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using Arblib
    using CGL
    using DataFrames
    using LaTeXStrings
    using OhMyThreads
    using Plots
    using Plots.PlotMeasures
    using PlutoUI

    setprecision(Arb, 128)

    ENV["LD_LIBRARY_PATH"] = "/home/joeldahne/Programs/capd/lib/"
end

# ╔═╡ b7903f77-8cd1-455f-aec7-be142ce71c4e
TableOfContents()

# ╔═╡ 6ea53a4c-7365-4599-b4e4-a4633c31063d
md"""
# Case I
"""

# ╔═╡ afe77be4-3366-45b5-b3eb-a1413c12c3d6
parameters_1 = [
    (j = 1, d = 1, ξ₁ = 10), # nothing gives default value
    (j = 2, d = 1, ξ₁ = 15),
    (j = 3, d = 1, ξ₁ = 20),
    (j = 4, d = 1, ξ₁ = 25),
    (j = 5, d = 1, ξ₁ = 40),
    (j = 6, d = 1, ξ₁ = 45),
    (j = 7, d = 1, ξ₁ = 60),
    (j = 8, d = 1, ξ₁ = 75),
]

# ╔═╡ a0609434-aeea-40f8-85e0-36ff26046aea
res_1 = OhMyThreads.tmap(parameters_1) do (j, d, ξ₁)
    μ, γ, κ, ϵ, ξ₁_actual, λ = CGL.sverak_params(Arb, j, d; ξ₁)
    CGL.G_solve_fix_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁_actual, λ), ξ₁_actual
end

# ╔═╡ f7d24e88-fe9a-4ff7-98dd-6db0d9d95736
critical_points_1 = OhMyThreads.tmap(res_1, parameters_1) do (sol, ξ₁), (j, d, _)
    μ, γ_real, γ_imag, κ = sol
    ϵ = zero(μ)
    λ = CGL.sverak_params(Arb, j, d)[6]
    CGL.count_critical_points(μ, Acb(γ_real, γ_imag), κ, ϵ, ξ₁, λ)
end

# ╔═╡ e2b0a9b6-31b3-443e-b1f2-872d0739b904
df_1 = let df = DataFrame()
    df.j = getfield.(parameters_1, :j)
    df.d = getfield.(parameters_1, :d)
    df.μ = getindex.(getindex.(res_1, 1), 1)
    df.γ_real = getindex.(getindex.(res_1, 1), 2)
    df.γ_imag = getindex.(getindex.(res_1, 1), 3)
    df.κ = getindex.(getindex.(res_1, 1), 4)
    df.ξ₁ = Float64.(getindex.(res_1, 2))
    df.num_critical_points =
        ifelse.(
            getindex.(critical_points_1, 1),
            length.(getindex.(critical_points_1, 2)),
            missing,
        )
    df
end

# ╔═╡ 285ba40d-98bd-4c23-9cd6-74fd3c5591df
md"""
# Case II
"""

# ╔═╡ 54c46505-cf47-404c-aa57-f331724951ee
parameters_2 = [
    (j = 1, d = 3, ξ₁ = 55),
    (j = 2, d = 3, ξ₁ = 140),
    (j = 3, d = 3, ξ₁ = nothing),
    (j = 4, d = 3, ξ₁ = nothing),
    (j = 5, d = 3, ξ₁ = nothing),
]

# ╔═╡ d836fdb5-8ea8-44a7-9c57-f4656545b27b
res_2 = OhMyThreads.tmap(parameters_2) do (j, d, ξ₁)
    μ, γ, κ, ϵ, ξ₁_actual, λ = CGL.sverak_params(Arb, j, d; ξ₁)
    CGL.G_solve_fix_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁_actual, λ), ξ₁_actual
end

# ╔═╡ cab3f05d-c949-490c-9367-6926732df899
critical_points_2 = OhMyThreads.tmap(res_2, parameters_2) do (sol, ξ₁), (j, d, _)
    μ, γ_real, γ_imag, κ = sol
    ϵ = zero(μ)
    λ = CGL.sverak_params(Arb, j, d)[6]
    CGL.count_critical_points(μ, Acb(γ_real, γ_imag), κ, ϵ, ξ₁, λ)
end

# ╔═╡ 2e75ba82-dbe8-472d-b7ab-e77175b1699a
df_2 = let df = DataFrame()
    df.j = getfield.(parameters_2, :j)
    df.d = getfield.(parameters_2, :d)
    df.μ = getindex.(getindex.(res_2, 1), 1)
    df.γ_real = getindex.(getindex.(res_2, 1), 2)
    df.γ_imag = getindex.(getindex.(res_2, 1), 3)
    df.κ = getindex.(getindex.(res_2, 1), 4)
    df.ξ₁ = Float64.(getindex.(res_2, 2))
    df.num_critical_points =
        ifelse.(
            getindex.(critical_points_2, 1),
            length.(getindex.(critical_points_2, 2)),
            missing,
        )
    df
end

# ╔═╡ 947a1fc3-6d56-4349-9cfa-1096a32f9241
md"""
# LaTeX output
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

# ╔═╡ 304cb8f0-b93c-438a-a327-5b78935244b5
md"""
# Temporary testing of solving
"""

# ╔═╡ 161ae97c-b0b5-4a4b-87fa-b5b80f9ba33c
# ╠═╡ disabled = true
#=╠═╡
res_tmp = let (j, d, ξ₁) = parameters_2[1]
    μ, γ, κ, ϵ, ξ₁_actual, λ = CGL.sverak_params(Arb, j, d; ξ₁)
    CGL.G_solve_fix_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁_actual, λ, verbose = true)
end
  ╠═╡ =#

# ╔═╡ 5b13d15d-d348-4f77-8472-123efc71836d
# ╠═╡ disabled = true
#=╠═╡
res_tmp = let (j, d, ξ₁) = parameters_2[1]
    μ, γ, κ, ϵ, ξ₁_actual, λ = CGL.sverak_params(Arb, j, d; ξ₁)
    CGL.G_solve_fix_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁_actual, λ, verbose = true)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═031b0a42-9f10-11ee-28d2-5957efedccab
# ╠═b7903f77-8cd1-455f-aec7-be142ce71c4e
# ╟─6ea53a4c-7365-4599-b4e4-a4633c31063d
# ╠═afe77be4-3366-45b5-b3eb-a1413c12c3d6
# ╠═a0609434-aeea-40f8-85e0-36ff26046aea
# ╠═f7d24e88-fe9a-4ff7-98dd-6db0d9d95736
# ╠═e2b0a9b6-31b3-443e-b1f2-872d0739b904
# ╟─285ba40d-98bd-4c23-9cd6-74fd3c5591df
# ╠═54c46505-cf47-404c-aa57-f331724951ee
# ╠═d836fdb5-8ea8-44a7-9c57-f4656545b27b
# ╠═cab3f05d-c949-490c-9367-6926732df899
# ╠═2e75ba82-dbe8-472d-b7ab-e77175b1699a
# ╟─947a1fc3-6d56-4349-9cfa-1096a32f9241
# ╟─eecc2ec2-ffb1-434a-a6a8-490672b56284
# ╟─a2449e4c-87f9-44ef-a1af-0eac9581f027
# ╠═36e2d4b8-60bd-418c-8311-328b1f6ed937
# ╠═e9d7b985-3e7c-4900-9d36-511d07ee6d6c
# ╟─304cb8f0-b93c-438a-a327-5b78935244b5
# ╠═161ae97c-b0b5-4a4b-87fa-b5b80f9ba33c
# ╠═5b13d15d-d348-4f77-8472-123efc71836d
