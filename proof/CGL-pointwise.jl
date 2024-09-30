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

# ╔═╡ 1628af7e-7e77-11ef-1ec8-55b8641b8458
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using Arblib
    using CGL
    using DataFrames
    using LaTeXStrings
    using OhMyThreads
    using Plots
    using PlutoUI

    setprecision(Arb, 128)
end

# ╔═╡ bcc838b0-078f-47c6-98be-c06e07ca1004
md"""
# Pointwise results for the CGL equation
This notebook presents the pointwise results for the CGL equation. It generates the figures for Section 9 in the paper. It primarily consists of precomputed data.
"""

# ╔═╡ 27579c86-563b-4844-b9ec-ed20f3fa8f92
TableOfContents()

# ╔═╡ 0e7a47c1-1c0b-49de-b99d-8a6531596852
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save CheckBox(default = false))
"""

# ╔═╡ 5c3d0478-1f01-4e18-bc8c-0e475578eee4
guidefontsize, tickfontsize = if save
    pgfplotsx()
    24, 24
else
    gr()
    11, 11
end

# ╔═╡ 7402af3a-23f7-4222-be67-456048c76b8d
md"""
## Load precomputed data
"""

# ╔═╡ e11054bb-e4aa-4420-9b3c-4e860793065a
branch_points_d1_epsilon = let
    dirname = "data/branch_points_d=1_fix_epsilon/"
    files = filter(endswith(".gz"), readdir(dirname))

    map(files) do file
        CGL.read_branch_points_csv(joinpath(dirname, file))
    end
end;

# ╔═╡ de0e538d-460e-49be-94e9-0c4fd6075d28
branch_points_d1_kappa = let
    dirname = "data/branch_points_d=1_fix_kappa/"
    files = filter(endswith(".gz"), readdir(dirname))

    map(files) do file
        CGL.read_branch_points_csv(joinpath(dirname, file))
    end
end;

# ╔═╡ da8565d6-b736-4541-a710-f002872cc847
branch_points_d3_epsilon = let
    dirname = "data/branch_points_d=3_fix_epsilon/"
    files = filter(endswith(".gz"), readdir(dirname))

    map(files) do file
        CGL.read_branch_points_csv(joinpath(dirname, file))
    end
end;

# ╔═╡ e768f98b-a1d4-4612-ae6d-9ada437e344f
branch_points_d3_kappa = let
    dirname = "data/branch_points_d=3_fix_kappa/"
    files = filter(endswith(".gz"), readdir(dirname))

    map(files) do file
        CGL.read_branch_points_csv(joinpath(dirname, file))
    end
end;

# ╔═╡ 83a33b9e-55df-48d7-9299-f6f4f75ca397
md"""
## Points per curve
"""

# ╔═╡ 5592fd67-e7ae-4c4d-87f5-2ecd1e4f44ac
nrow.(branch_points_d1_epsilon)

# ╔═╡ fe9fadd6-bb37-4691-962b-07fc61b12739
nrow.(branch_points_d1_kappa)

# ╔═╡ fd1d64c9-5c70-4e08-8541-95d25606fad9
nrow.(branch_points_d3_epsilon)

# ╔═╡ 5f7c7f39-2e18-434c-9017-85245e352508
nrow.(branch_points_d3_kappa)

# ╔═╡ 7c641571-e368-42c9-bf29-6754c3a41071
md"""
## Plot data
"""

# ╔═╡ 9ae41083-81c2-4607-8209-f353bee948d0
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa";
        guidefontsize,
        tickfontsize,
    )
    foreach(branch_points_d1_epsilon) do data
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
    save && savefig(pl, "figures/CGL-pointwise-d=1-epsilon.pdf")
    pl
end

# ╔═╡ dccce64d-3184-433a-8ad2-ee45876bfee4
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa";
        guidefontsize,
        tickfontsize,
    )
    foreach(branch_points_d1_kappa) do data
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
    save && savefig(pl, "figures/CGL-pointwise-d=1-kappa.pdf")
    pl
end

# ╔═╡ 88f14090-9885-4e1b-95b9-ebe3c8142eea
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa";
        guidefontsize,
        tickfontsize,
    )
    foreach(branch_points_d3_epsilon) do data
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
    save && savefig(pl, "figures/CGL-pointwise-d=3-epsilon.pdf")
    pl
end

# ╔═╡ 449ba275-2c0d-4629-bb0d-50fb215ce0ac
let pl = plot(
        legend = :none,
        xlabel = L"\epsilon",
        ylabel = L"\kappa";
        guidefontsize,
        tickfontsize,
    )
    foreach(branch_points_d3_kappa) do data
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
    save && savefig(pl, "figures/CGL-pointwise-d=3-kappa.pdf")
    pl
end

# ╔═╡ 19d93cd7-d9e4-495d-83e6-0f5b98138a7d
let pl = plot(
        legend = :none,
        xlabel = L"\kappa",
        ylabel = L"\xi_1";
        guidefontsize,
        tickfontsize,
    )

    scatter!(
        pl,
        Float64.(branch_points_d3_kappa[5].κ),
        Float64.(branch_points_d3_kappa[5].ξ₁),
        markerstrokewidth = 0,
    )

    pl
end

# ╔═╡ Cell order:
# ╟─bcc838b0-078f-47c6-98be-c06e07ca1004
# ╠═1628af7e-7e77-11ef-1ec8-55b8641b8458
# ╠═27579c86-563b-4844-b9ec-ed20f3fa8f92
# ╟─0e7a47c1-1c0b-49de-b99d-8a6531596852
# ╟─5c3d0478-1f01-4e18-bc8c-0e475578eee4
# ╟─7402af3a-23f7-4222-be67-456048c76b8d
# ╠═e11054bb-e4aa-4420-9b3c-4e860793065a
# ╠═de0e538d-460e-49be-94e9-0c4fd6075d28
# ╠═da8565d6-b736-4541-a710-f002872cc847
# ╠═e768f98b-a1d4-4612-ae6d-9ada437e344f
# ╟─83a33b9e-55df-48d7-9299-f6f4f75ca397
# ╠═5592fd67-e7ae-4c4d-87f5-2ecd1e4f44ac
# ╠═fe9fadd6-bb37-4691-962b-07fc61b12739
# ╠═fd1d64c9-5c70-4e08-8541-95d25606fad9
# ╠═5f7c7f39-2e18-434c-9017-85245e352508
# ╟─7c641571-e368-42c9-bf29-6754c3a41071
# ╟─9ae41083-81c2-4607-8209-f353bee948d0
# ╟─dccce64d-3184-433a-8ad2-ee45876bfee4
# ╟─88f14090-9885-4e1b-95b9-ebe3c8142eea
# ╟─449ba275-2c0d-4629-bb0d-50fb215ce0ac
# ╠═19d93cd7-d9e4-495d-83e6-0f5b98138a7d
