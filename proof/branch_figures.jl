### A Pluto.jl notebook ###
# v0.19.42

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

# ╔═╡ affd856e-dd45-11ee-1936-75eeb3c8c5bd
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using Arblib
    using CGL
    using DataFrames
    using LaTeXStrings
    using OhMyThreads: tmap
    using Plots
    using Plots.PlotMeasures
    using PlutoUI
    using StaticArrays
    using IntervalArithmetic

    import CGL: G, G_jacobian_kappa
    setprecision(Arb, 128)
    pgfplotsx()
end

# ╔═╡ 6fa482f8-ac70-4e81-83e7-57906ce13c79
function IntervalArithmetic.interval(x::Arb)
    if isnan(Arblib.midref(x))
        return nai(Float64)
    else
        return interval(Float64, getinterval(BigFloat, x)...)
    end
end

# ╔═╡ f447b12c-943b-45fa-872a-3ca9b81b281f
function IntervalArithmetic.interval((x, y)::NTuple{2,Arb})
    if isnan(Arblib.midref(x)) || isnan(Arblib.midref(y))
        return nai(Float64)
    else
        return interval(Float64, BigFloat(lbound(x)), BigFloat(ubound(y)))
    end
end

# ╔═╡ 9898b8dd-f08d-4c0e-95ef-86710e9d7d51
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save CheckBox(default = false))
"""

# ╔═╡ 5875edd0-fbe5-47a5-86b3-f5226e9d3e7d
md"""
## Compute branches
"""

# ╔═╡ cee4af3d-1f29-4021-9546-b2efa3e90d1c
branches_d1 = tmap(1:8, chunksize = 1) do j
    μ, κ, ϵ, ω, λ = CGL.CGLBranch.sverak_initial(j, 1)
    CGL.CGLBranch.branch_epsilon(μ, κ, ϵ, ω, λ)
end

# ╔═╡ 1825a189-7a6f-437c-bab7-23b822c30241
branches_d3 = tmap(1:5, chunksize = 1) do j
    μ, κ, ϵ, ω, λ = CGL.CGLBranch.sverak_initial(j, 3)
    CGL.CGLBranch.branch_epsilon(μ, κ, ϵ, ω, λ)
end

# ╔═╡ 1902e6a7-aa08-43dc-97d2-a985c144e2ad
md"""
## Load precomputed data
"""

# ╔═╡ 637c3a83-e7b7-4187-9172-9fe748e398fe
branches_points_epsilon_d1 = let
    dirname = "data/branches_points_epsilon_d=1/"
    files = readdir(dirname)

    map(files) do file
        CGL.read_branch_points_csv(joinpath(dirname, file))
    end
end;

# ╔═╡ be6126e7-cab3-4733-8adc-f5f32fe90e44
branches_points_kappa_d1 = let
    dirname = "data/branches_points_kappa_d=1/"
    files = readdir(dirname)

    map(files) do file
        CGL.read_branch_points_csv(joinpath(dirname, file))
    end
end;

# ╔═╡ d6139ee8-fb02-4b55-becb-b7b92df5bca7
branches_points_epsilon_d3 = let
    dirname = "data/branches_points_epsilon_d=3/"
    files = readdir(dirname)

    map(files) do file
        CGL.read_branch_points_csv(joinpath(dirname, file))
    end
end;

# ╔═╡ 94408901-278c-4a32-9543-7538e8490f77
branches_points_kappa_d3 = let
    dirname = "data/branches_points_kappa_d=3/"
    files = readdir(dirname)

    map(files) do file
        CGL.read_branch_points_csv(joinpath(dirname, file))
    end
end;

# ╔═╡ 97712600-eea3-48a4-afc2-092af4e746d5
parameters, data_top, data_turn, data_bottom, data_connection_points =
    CGL.read_proof_witness("data/branch_j=3_d=1/");

# ╔═╡ 59277a78-9050-41ce-9d8b-af65b6aa883b
all_branches = [CGL.read_proof_witness("data/branch_j=$(j)_d=1/") for j = 1:8];

# ╔═╡ e6be97f8-1dad-455b-898e-d97c8697a3d5
md"""
## Plot branches
"""

# ╔═╡ 059599ed-4b68-4b32-83e2-8d9716f206f6
md"""
### Numerical approximation
"""

# ╔═╡ 876da690-ef1d-4f34-ae5b-5dd3bdef89d2
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        guidefontsize = 24,
        tickfontsize = 24,
    )
    foreach(br -> plot!(pl, br.param, br.κ), branches_d1)
    save && savefig(pl, "figures/branches-d=1-numerical.pdf")
    pl
end

# ╔═╡ 25af589a-6ca1-47ec-b933-dd2addc4eb3b
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        guidefontsize = 24,
        tickfontsize = 24,
    )
    foreach(br -> plot!(pl, br.param, br.κ), branches_d3)
    save && savefig(pl, "figures/branches-d=3-numerical.pdf")
    pl
end

# ╔═╡ 875f10b0-3280-4a1e-8f24-51330ee2a7f0
md"""
### Pointwise verification
"""

# ╔═╡ 37a33989-23a2-45c1-856b-7bd181f69da5
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        guidefontsize = 24,
        tickfontsize = 24,
    )
    foreach(branches_points_epsilon_d1) do data
        is_finite = isfinite.(data.μ)
        scatter!(
            pl,
            Float64.(data.ϵ[is_finite]),
            Float64.(data.κ[is_finite]),
            color = :blue,
            markerstrokecolor = :blue,
            markersize = 1,
        )
        scatter!(
            pl,
            Float64.(data.ϵ[.!is_finite]),
            Float64.(data.κ₀[.!is_finite]),
            color = :red,
            markerstrokecolor = :red,
            markersize = 1,
        )
    end
    save && savefig(pl, "figures/branches-points-epsilon-d=1.pdf")
    pl
end

# ╔═╡ 26810775-aa38-4f51-8b76-9209fc742371
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        guidefontsize = 24,
        tickfontsize = 24,
    )
    foreach(branches_points_kappa_d1) do data
        is_finite = isfinite.(data.μ)
        scatter!(
            pl,
            Float64.(data.ϵ[is_finite]),
            Float64.(data.κ[is_finite]),
            color = :blue,
            markerstrokecolor = :blue,
            markersize = 1,
        )
        scatter!(
            pl,
            Float64.(data.ϵ₀[.!is_finite]),
            Float64.(data.κ[.!is_finite]),
            color = :red,
            markerstrokecolor = :red,
            markersize = 1,
        )
    end
    save && savefig(pl, "figures/branches-points-kappa-d=1.pdf")
    pl
end

# ╔═╡ ad97d3df-e931-4a09-ad4c-7545bf9da7a1
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        guidefontsize = 24,
        tickfontsize = 24,
    )
    foreach(branches_points_epsilon_d3) do data
        is_finite = isfinite.(data.μ)
        scatter!(
            pl,
            Float64.(data.ϵ[is_finite]),
            Float64.(data.κ[is_finite]),
            color = :blue,
            markerstrokecolor = :blue,
            markersize = 1,
        )
        scatter!(
            pl,
            Float64.(data.ϵ[.!is_finite]),
            Float64.(data.κ₀[.!is_finite]),
            color = :red,
            markerstrokecolor = :red,
            markersize = 1,
        )
    end
    save && savefig(pl, "figures/branches-points-epsilon-d=3.pdf")
    pl
end

# ╔═╡ 70570c41-ab75-46c9-8780-6ee7ad1772a1
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        guidefontsize = 24,
        tickfontsize = 24,
    )
    foreach(branches_points_kappa_d3) do data
        is_finite = isfinite.(data.μ)
        scatter!(
            pl,
            Float64.(data.ϵ[is_finite]),
            Float64.(data.κ[is_finite]),
            color = :blue,
            markerstrokecolor = :blue,
            markersize = 1,
        )
        scatter!(
            pl,
            Float64.(data.ϵ₀[.!is_finite]),
            Float64.(data.κ[.!is_finite]),
            color = :red,
            markerstrokecolor = :red,
            markersize = 1,
        )
    end
    save && savefig(pl, "figures/branches-points-kappa-d=3.pdf")
    pl
end

# ╔═╡ fbe496b0-0e4c-4113-b1a1-c02fe05d1ae0
md"""
## Continuous verification
"""

# ╔═╡ b575f9f8-32f4-414c-8b14-6c2d245b07c4
md"""
### $d = 1$, $j = 3$
"""

# ╔═╡ 124fcb05-aa32-4149-87d7-8755d58de880
let
    pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        guidefontsize = 24,
        tickfontsize = 24,
    )

    top_indices = round.(Int, range(1, nrow(data_top), length = 2000))
    turn_indices = round.(Int, range(1, nrow(data_turn), length = 2000))
    bottom_indices = round.(Int, range(1, nrow(data_bottom), length = 2000))

    plot!(
        pl,
        Float64.(data_top.ϵ_lower[top_indices]),
        Float64.(data_top.κ_exists[top_indices]),
        color = :green,
        label = "",
    )
    plot!(
        pl,
        Float64.(data_turn.ϵ_exists[turn_indices]),
        Float64.(data_turn.κ_lower[turn_indices]),
        color = :blue,
        label = "",
    )
    plot!(
        pl,
        Float64.(data_bottom.ϵ_lower[bottom_indices]),
        Float64.(data_bottom.κ_exists[bottom_indices]),
        color = :green,
        label = "",
    )

    scatter!(
        Float64.(data_connection_points.ϵ),
        Float64.(data_connection_points.κ_exists),
        markersize = 3,
        color = :red,
        label = "Connection points",
    )
    save && savefig(pl, "figures/branch-j=3-d=1.pdf")
    pl
end

# ╔═╡ 1be510f4-f5d1-43fc-a494-e4b9fee552e8
let data = data_top[1:50, :]
    ϵs = tuple.(data.ϵ_lower, data.ϵ_upper)
    κs_uniq = data.κ_uniq
    κs_exists = data.κ_exists

    xticks = [0.0, 1e-4, 2e-4]
    yticks = [0.34662, 0.34675]

    pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        xticks = xticks,
        yticks = (yticks, string.(yticks)),
        guidefontsize = 24,
        tickfontsize = 24,
    )

    plot!(pl, vcat.(interval.(ϵs), interval.(κs_uniq)), linecolor = :red, color = :red)
    plot!(
        pl,
        vcat.(interval.(ϵs), interval.(κs_exists)),
        linecolor = :green,
        color = :green,
    )
    save && savefig(pl, "figures/branch-j=3-d=1-start.pdf")
    pl
end

# ╔═╡ 206b74fa-60e9-4b17-b411-2b0db1a5e641
let
    point = SVector(
        data_connection_points.μ_exists[1],
        real(data_connection_points.γ_exists[1]),
        imag(data_connection_points.γ_exists[1]),
        data_connection_points.κ_exists[1],
        data_connection_points.ϵ[1],
    )

    i = findfirst(1:nrow(data_turn)) do i
        all(
            Arblib.contains_interior.(
                SVector(
                    data_turn.μ_uniq[i],
                    real(data_turn.γ_uniq[i]),
                    imag(data_turn.γ_uniq[i]),
                    Arb((data_turn.κ_lower[i], data_turn.κ_upper[i])),
                    data_turn.ϵ_uniq[i],
                ),
                point,
            ),
        )
    end

    ϵs_top = tuple.(data_top.ϵ_lower[end-3:end], data_top.ϵ_upper[end-3:end])
    κs_uniq_top = data_top.κ_uniq[end-3:end]
    κs_exists_top = data_top.κ_exists[end-3:end]

    ϵs_uniq_turn = data_turn.ϵ_uniq[i-3:i+3]
    ϵs_exists_turn = data_turn.ϵ_exists[i-3:i+3]
    κs_turn = tuple.(data_turn.κ_lower[i-3:i+3], data_turn.κ_upper[i-3:i+3])

    xticks = [0.022786, 0.02279]
    xlims_diff = 0.1(xticks[2] - xticks[1])
    xlims = xticks .+ xlims_diff * [-1, 1]
    yticks = [0.327434, 0.32745]
    ylims_diff = 0.1(yticks[2] - yticks[1])
    ylims = yticks .+ ylims_diff * [-1, 1]

    pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        xlims = xlims,
        ylims = ylims,
        xticks = (xticks, string.(xticks)),
        yticks = (yticks, string.(yticks)),
        guidefontsize = 24,
        tickfontsize = 24,
    )

    plot!(
        pl,
        vcat.(interval.(ϵs_exists_turn), interval.(κs_turn)),
        linecolor = :blue,
        color = :blue,
    )

    plot!(
        pl,
        vcat.(interval.(ϵs_top), interval.(κs_exists_top)),
        linecolor = :green,
        color = :green,
    )

    scatter!(
        pl,
        Float64.(data_connection_points.ϵ[1:1]),
        Float64.(data_connection_points.κ_exists[1:1]),
        markersize = 3,
        color = :red,
    )
    save && savefig(pl, "figures/branch-j=3-d=1-connection-top.pdf")
    pl
end

# ╔═╡ 66901c61-0551-42a0-a6b5-cc6d61cf4f51
let
    point = SVector(
        data_connection_points.μ_exists[2],
        real(data_connection_points.γ_exists[2]),
        imag(data_connection_points.γ_exists[2]),
        data_connection_points.κ_exists[2],
        data_connection_points.ϵ[2],
    )

    i = findlast(1:nrow(data_turn)) do i
        all(
            Arblib.contains_interior.(
                SVector(
                    data_turn.μ_uniq[i],
                    real(data_turn.γ_uniq[i]),
                    imag(data_turn.γ_uniq[i]),
                    Arb((data_turn.κ_lower[i], data_turn.κ_upper[i])),
                    data_turn.ϵ_uniq[i],
                ),
                point,
            ),
        )
    end

    ϵs_bottom =
        interval.(Float64, tuple.(data_bottom.ϵ_lower[1:10], data_bottom.ϵ_upper[1:10]))
    κs_exists_bottom = interval.(data_bottom.κ_exists[1:10])

    ϵs_exists_turn = interval.(data_turn.ϵ_exists[i-5:i+5])
    κs_turn =
        interval.(Float64, tuple.(data_turn.κ_lower[i-5:i+5], data_turn.κ_upper[i-5:i+5]))

    xticks = [0.0324785, 0.03248]
    xlims_diff = 0.1(xticks[2] - xticks[1])
    xlims = xticks .+ xlims_diff * [-1, 1]
    yticks = [0.19438, 0.19441]
    ylims_diff = 0.1(yticks[2] - yticks[1])
    ylims = yticks .+ ylims_diff * [-1, 1]

    pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        xlims = xlims,
        ylims = ylims,
        xticks = (xticks, string.(xticks)),
        yticks = (yticks, string.(yticks)),
        guidefontsize = 24,
        tickfontsize = 24,
    )

    plot!(pl, vcat.(ϵs_exists_turn, κs_turn), linecolor = :blue, color = :blue)

    plot!(pl, vcat.(ϵs_bottom, κs_exists_bottom), linecolor = :green, color = :green)

    scatter!(
        pl,
        Float64.(data_connection_points.ϵ[2:2]),
        Float64.(data_connection_points.κ_exists[2:2]),
        markersize = 3,
        color = :red,
    )
    save && savefig(pl, "figures/branch-j=3-d=1-connection-bottom.pdf")
    pl
end

# ╔═╡ 8ea15d67-da43-46ad-86a6-383a4e587abf
md"""
### $d = 1$, all $j$
"""

# ╔═╡ 1dd55d16-e27d-4bc9-8f01-fcd321a5b58d
let
    pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa",
        guidefontsize = 24,
        tickfontsize = 24,
    )

    for (parameters, data_top, data_turn, data_bottom, data_connection_points) in
        all_branches

        top_indices = round.(Int, range(1, nrow(data_top), length = 1000))
        plot!(
            pl,
            Float64.(data_top.ϵ_lower[top_indices]),
            Float64.(data_top.κ_exists[top_indices]),
            color = :green,
            label = "",
        )

        if !isnothing(data_turn)
            turn_indices = round.(Int, range(1, nrow(data_turn), length = 1000))
            plot!(
                pl,
                Float64.(data_turn.ϵ_exists[turn_indices]),
                Float64.(data_turn.κ_lower[turn_indices]),
                color = :blue,
                label = "",
            )
        end

        if !isnothing(data_bottom)
            bottom_indices = round.(Int, range(1, nrow(data_bottom), length = 1000))
            plot!(
                pl,
                Float64.(data_bottom.ϵ_lower[bottom_indices]),
                Float64.(data_bottom.κ_exists[bottom_indices]),
                color = :green,
                label = "",
            )
        end

        scatter!(
            Float64.(data_connection_points.ϵ),
            Float64.(data_connection_points.κ_exists),
            markersize = 3,
            color = :red,
            label = "Connection points",
        )
    end
    save && savefig(pl, "figures/branches-d=1.pdf")
    pl
end

# ╔═╡ e36b9e3f-6f92-43a5-a830-a5a8049090a1
md"""
## Failure of continuation
"""

# ╔═╡ 88fe07f2-0e7a-4c74-9603-4c491788762a
let pl = plot(legend = :none, axis = ([], false), grid = false)
    plot!(pl, range(-0.1, 0.1, 1000), x -> 0.998cos(x), color = :blue)
    plot!(pl, range(-0.1, 0.1, 1000), x -> 2 - cos(x), color = :blue)
    plot!(pl, [interval(-0.01, 0), interval(0.9965, 0.9995)], color = :yellow)
    plot!(pl, [interval(-0.01, 0), interval(0.9968, 0.9992)], color = :green)
    plot!(pl, [interval(0, 0.01), interval(0.9985, 1.0015)], color = :yellow)
    plot!(pl, [interval(0, 0.01), interval(0.9988, 1.0012)], color = :green)
    pl
end

# ╔═╡ Cell order:
# ╠═affd856e-dd45-11ee-1936-75eeb3c8c5bd
# ╠═6fa482f8-ac70-4e81-83e7-57906ce13c79
# ╠═f447b12c-943b-45fa-872a-3ca9b81b281f
# ╟─9898b8dd-f08d-4c0e-95ef-86710e9d7d51
# ╟─5875edd0-fbe5-47a5-86b3-f5226e9d3e7d
# ╠═cee4af3d-1f29-4021-9546-b2efa3e90d1c
# ╠═1825a189-7a6f-437c-bab7-23b822c30241
# ╟─1902e6a7-aa08-43dc-97d2-a985c144e2ad
# ╠═637c3a83-e7b7-4187-9172-9fe748e398fe
# ╠═be6126e7-cab3-4733-8adc-f5f32fe90e44
# ╠═d6139ee8-fb02-4b55-becb-b7b92df5bca7
# ╠═94408901-278c-4a32-9543-7538e8490f77
# ╠═97712600-eea3-48a4-afc2-092af4e746d5
# ╠═59277a78-9050-41ce-9d8b-af65b6aa883b
# ╟─e6be97f8-1dad-455b-898e-d97c8697a3d5
# ╟─059599ed-4b68-4b32-83e2-8d9716f206f6
# ╟─876da690-ef1d-4f34-ae5b-5dd3bdef89d2
# ╟─25af589a-6ca1-47ec-b933-dd2addc4eb3b
# ╟─875f10b0-3280-4a1e-8f24-51330ee2a7f0
# ╟─37a33989-23a2-45c1-856b-7bd181f69da5
# ╟─26810775-aa38-4f51-8b76-9209fc742371
# ╟─ad97d3df-e931-4a09-ad4c-7545bf9da7a1
# ╟─70570c41-ab75-46c9-8780-6ee7ad1772a1
# ╟─fbe496b0-0e4c-4113-b1a1-c02fe05d1ae0
# ╟─b575f9f8-32f4-414c-8b14-6c2d245b07c4
# ╟─124fcb05-aa32-4149-87d7-8755d58de880
# ╟─1be510f4-f5d1-43fc-a494-e4b9fee552e8
# ╟─206b74fa-60e9-4b17-b411-2b0db1a5e641
# ╟─66901c61-0551-42a0-a6b5-cc6d61cf4f51
# ╟─8ea15d67-da43-46ad-86a6-383a4e587abf
# ╟─1dd55d16-e27d-4bc9-8f01-fcd321a5b58d
# ╟─e36b9e3f-6f92-43a5-a830-a5a8049090a1
# ╠═88fe07f2-0e7a-4c74-9603-4c491788762a
