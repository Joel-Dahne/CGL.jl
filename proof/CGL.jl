### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(
                Base.UUID("6e696c72-6542-2067-7265-42206c756150"),
                "AbstractPlutoDingetjes",
            )].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 02cab23a-7e77-11ef-2d7a-73b1a4cefdc3
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
    using StatsPlots

    setprecision(Arb, 128)

    nothing
end

# ╔═╡ 2541f0ac-36df-40e3-acd3-0bb33412875e
md"""
# Results for the CGL equation
This notebook presents the results for the CGL equation. It generates the figures and tables for Section 4 and 6.2 in the paper. It primarily conists of precomputed data.
"""

# ╔═╡ 165c76cc-f6d1-4f51-b159-bb1e489fccd4
TableOfContents()

# ╔═╡ b2020ad6-911d-4c93-901c-56ad4328deea
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save CheckBox(default = false))
"""

# ╔═╡ 8769a571-1479-4baf-bf66-89e7d50257b4
guidefontsize, tickfontsize = if save
    pgfplotsx()
    24, 24
else
    gr()
    11, 11
end

# ╔═╡ c22c2eaa-86ea-49f0-8f65-0a315f09e01a
md"""
## Compute numerical approximations of branches
"""

# ╔═╡ a5d2fdab-434a-4a60-a16f-5bf4b99c6fbe
branches_approximation_d1 = tmap(1:8, chunksize = 1) do j
    CGL.CGLBranch.branch_epsilon(CGL.CGLBranch.sverak_initial(j, 1)...)
end

# ╔═╡ 75595e55-06d7-435d-9c8b-3f078698a8f1
branches_approximation_d3 = tmap(1:5, chunksize = 1) do j
    CGL.CGLBranch.branch_epsilon(CGL.CGLBranch.sverak_initial(j, 3)...)
end

# ╔═╡ 77224e16-7066-47da-935a-a74df7da4559
md"""
## Read precomputed data
"""

# ╔═╡ ca542d21-914a-43c8-a436-cf18400dda52
branches_d1 = [CGL.read_proof_witness("data/branch_d=1_j=$(j)/") for j = 1:8];

# ╔═╡ 7ec01aaa-08d5-4f29-ba9e-5cde5a3551e6
branches_d3 = [CGL.read_proof_witness("data/branch_d=3_j=$(j)/") for j = 1:4];

# ╔═╡ 816cb3d9-bf33-40d6-9a9f-d24ac8afb596
md"""
## Check precomputed data
"""

# ╔═╡ 1b7bc518-87fe-420a-982f-920650d3772c
for (j, (parameters, data_top, data_turn, data_bottom, data_connection_points)) in
    enumerate(branches_d1)
    @info "Checking j = $j"
    CGL.check_proof_witness(
        parameters,
        data_top,
        data_turn,
        data_bottom,
        data_connection_points,
        assert_critical_points = true,
    )
end

# ╔═╡ 589afeb9-ed71-4778-82f3-b2b1365fb0e8
for (j, (parameters, data_top, data_turn, data_bottom, data_connection_points)) in
    enumerate(branches_d3)
    @info "Checking j = $j"
    CGL.check_proof_witness(
        parameters,
        data_top,
        data_turn,
        data_bottom,
        data_connection_points,
    )
end

# ╔═╡ b72016ec-e78e-4b14-bfcc-5bfbec1d5771
md"""
## Plot data
"""

# ╔═╡ 9b7d37a0-a45b-4c29-beb1-f791d66b078a
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa";
        guidefontsize,
        tickfontsize,
    )
    foreach(br -> plot!(pl, br.param, br.κ, linestyle = :dot), branches_approximation_d1)
    save && savefig(pl, "figures/CGL-approximation-d=1.pdf")
    pl
end

# ╔═╡ feb2eb52-30bd-4b37-99dc-5884c0212c1a
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa";
        guidefontsize,
        tickfontsize,
    )
    foreach(br -> plot!(pl, br.param, br.κ, linestyle = :dot), branches_approximation_d3)
    save && savefig(pl, "figures/CGL-approximation-d=3.pdf")
    pl
end

# ╔═╡ 9b51a77e-92b7-4324-8f56-158a62367087
let
    pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa";
        guidefontsize,
        tickfontsize,
    )

    for (j, (parameters, data_top, data_turn, data_bottom, data_connection_points)) in
        enumerate(branches_d1)

        top_indices = round.(Int, range(1, nrow(data_top), length = 1000))
        plot!(
            pl,
            Float64.(data_top.ϵ_lower[top_indices]),
            Float64.(data_top.κ_exists[top_indices]),
            color = j,
            label = "",
        )

        if !isnothing(data_turn)
            turn_indices = round.(Int, range(1, nrow(data_turn), length = 1000))
            plot!(
                pl,
                Float64.(data_turn.ϵ_exists[turn_indices]),
                Float64.(data_turn.κ_lower[turn_indices]),
                color = j,
                label = "",
            )
        end

        if !isnothing(data_bottom)
            bottom_indices = round.(Int, range(1, nrow(data_bottom), length = 1000))
            plot!(
                pl,
                Float64.(data_bottom.ϵ_lower[bottom_indices]),
                Float64.(data_bottom.κ_exists[bottom_indices]),
                color = j,
                label = "",
            )
        end
    end
    save && savefig(pl, "figures/CGL-d=1.pdf")
    pl
end

# ╔═╡ a5eb6957-aea3-4a49-889e-027cc9fcd611
let
    pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa";
        guidefontsize,
        tickfontsize,
    )

    foreach(branches_approximation_d3[eachindex(branches_d3)]) do br
        plot!(pl, br.param, br.κ, linestyle = :dot)
    end

    for (j, (parameters, data_top, data_turn, data_bottom, data_connection_points)) in
        enumerate(branches_d3)

        if !isnothing(data_top)
            top_indices = round.(Int, range(1, nrow(data_top), length = 1000))
            plot!(
                pl,
                Float64.(data_top.ϵ_lower[top_indices]),
                Float64.(data_top.κ_exists[top_indices]),
                color = j,
                label = "",
                linewidth = 2,
            )
        end

        if !isnothing(data_turn)
            turn_indices = round.(Int, range(1, nrow(data_turn), length = 1000))
            plot!(
                pl,
                Float64.(data_turn.ϵ_exists[turn_indices]),
                Float64.(data_turn.κ_lower[turn_indices]),
                color = j,
                label = "",
                linewidth = 2,
            )
        end

        if !isnothing(data_bottom)
            bottom_indices = round.(Int, range(1, nrow(data_bottom), length = 1000))
            plot!(
                pl,
                Float64.(data_bottom.ϵ_lower[bottom_indices]),
                Float64.(data_bottom.κ_exists[bottom_indices]),
                color = j,
                label = "",
                linewidth = 2,
            )
        end
    end
    save && savefig(pl, "figures/CGL-d=3.pdf")
    pl
end

# ╔═╡ f717a6ba-063d-4e11-b504-ae0a4ec150fd
md"""
## Computational time
"""

# ╔═╡ 4729a916-81f6-4b01-95b2-fa2c7a81c30e
cores = 256 # For converting from seconds to core hours

# ╔═╡ 545a05dc-d87e-4a36-bddc-b0b1e0ef204d
let branches = branches_d1
    runtimes_top =
        cores / 3600 * map(branches) do branch
            branch[1].runtime_existence_top + branch[1].runtime_continuation_top
        end
    runtimes_turn =
        cores / 3600 * map(branches) do branch
            branch[1].runtime_existence_turn + branch[1].runtime_continuation_turn
        end
    runtimes_bottom =
        cores / 3600 * map(branches) do branch
            branch[1].runtime_existence_bottom + branch[1].runtime_continuation_bottom
        end

    pl = groupedbar(
        hcat(runtimes_top, runtimes_turn, runtimes_bottom),
        bar_position = :stack,
        labels = ["Top" "Turn" "Bottom"],
        xticks = eachindex(branches),
        xlabel = L"j",
        ylabel = "Core hours",
        legend = :topright;
        guidefontsize,
        tickfontsize,
    )

    save && savefig(pl, "figures/CGL-runtime-d=1.pdf")

    pl
end

# ╔═╡ 06883f44-7744-4a02-a301-4a9a327127e7
let branches = branches_d1
    runtimes_critica_points_top = cores / 3600 * map(branches) do branch
        branch[1].runtime_critical_points_top
    end
    runtimes_critica_points_turn = cores / 3600 * map(branches) do branch
        branch[1].runtime_critical_points_turn
    end
    runtimes_critica_points_bottom =
        cores / 3600 * map(branches) do branch
            branch[1].runtime_critical_points_bottom
        end

    pl = groupedbar(
        hcat(
            runtimes_critica_points_top,
            runtimes_critica_points_turn,
            runtimes_critica_points_bottom,
        ),
        bar_position = :stack,
        labels = ["Top" "Turn" "Bottom"],
        xticks = eachindex(branches),
        xlabel = L"j",
        ylabel = "Core hours",
        legend = :topright;
        guidefontsize,
        tickfontsize,
    )

    save && savefig(pl, "figures/CGL-runtime-critical-points-d=1.pdf")

    pl
end

# ╔═╡ 9c02cabe-9dce-470a-be84-3f31e4a7a21e
let stats = DataFrame()
    runtimes_top =
        cores / 3600 * map(branches_d3) do branch
            branch[1].runtime_existence_top + branch[1].runtime_continuation_top
        end
    runtimes_turn =
        cores / 3600 * map(branches_d3) do branch
            branch[1].runtime_existence_turn + branch[1].runtime_continuation_turn
        end
    runtimes_bottom =
        cores / 3600 * map(branches_d3) do branch
            branch[1].runtime_existence_bottom + branch[1].runtime_continuation_bottom
        end

    pl = groupedbar(
        hcat(runtimes_top, runtimes_turn, runtimes_bottom),
        bar_position = :stack,
        labels = ["Top" "Turn" "Bottom"],
        xticks = eachindex(branches_d3),
        xlabel = L"j",
        ylabel = "Core hours",
        legend = :topright;
        guidefontsize,
        tickfontsize,
    )

    save && savefig(pl, "figures/CGL-runtime-d=3.pdf")

    pl
end

# ╔═╡ Cell order:
# ╟─2541f0ac-36df-40e3-acd3-0bb33412875e
# ╠═02cab23a-7e77-11ef-2d7a-73b1a4cefdc3
# ╠═165c76cc-f6d1-4f51-b159-bb1e489fccd4
# ╟─b2020ad6-911d-4c93-901c-56ad4328deea
# ╟─8769a571-1479-4baf-bf66-89e7d50257b4
# ╟─c22c2eaa-86ea-49f0-8f65-0a315f09e01a
# ╠═a5d2fdab-434a-4a60-a16f-5bf4b99c6fbe
# ╠═75595e55-06d7-435d-9c8b-3f078698a8f1
# ╟─77224e16-7066-47da-935a-a74df7da4559
# ╠═ca542d21-914a-43c8-a436-cf18400dda52
# ╠═7ec01aaa-08d5-4f29-ba9e-5cde5a3551e6
# ╟─816cb3d9-bf33-40d6-9a9f-d24ac8afb596
# ╠═1b7bc518-87fe-420a-982f-920650d3772c
# ╠═589afeb9-ed71-4778-82f3-b2b1365fb0e8
# ╟─b72016ec-e78e-4b14-bfcc-5bfbec1d5771
# ╟─9b7d37a0-a45b-4c29-beb1-f791d66b078a
# ╟─feb2eb52-30bd-4b37-99dc-5884c0212c1a
# ╟─9b51a77e-92b7-4324-8f56-158a62367087
# ╟─a5eb6957-aea3-4a49-889e-027cc9fcd611
# ╟─f717a6ba-063d-4e11-b504-ae0a4ec150fd
# ╠═4729a916-81f6-4b01-95b2-fa2c7a81c30e
# ╟─545a05dc-d87e-4a36-bddc-b0b1e0ef204d
# ╟─06883f44-7744-4a02-a301-4a9a327127e7
# ╟─9c02cabe-9dce-470a-be84-3f31e4a7a21e
