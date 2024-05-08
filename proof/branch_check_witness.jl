### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 2f550824-d714-11ee-1c46-9b93dd3db181
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using Arblib
    using CGL
    using DataFrames
    using LaTeXStrings
    using Plots
    using Plots.PlotMeasures
    using PlutoUI
    using StaticArrays
    using IntervalArithmetic

    import CGL: G, G_jacobian_kappa
    setprecision(Arb, 128)
end

# ╔═╡ ec89485a-4c4d-45d3-b2fe-1a1a1d61f147
function IntervalArithmetic.interval(x::Arb)
    if isnan(Arblib.midref(x))
        return nai(Float64)
    else
        return interval(Float64, getinterval(BigFloat, x)...)
    end
end

# ╔═╡ b073c9ef-0dcd-4496-b72b-f36fac091474
function IntervalArithmetic.interval((x, y)::NTuple{2,Arb})
    if isnan(Arblib.midref(x)) || isnan(Arblib.midref(y))
        return nai(Float64)
    else
        return interval(Float64, BigFloat(lbound(x)), BigFloat(ubound(y)))
    end
end

# ╔═╡ ed29ab9a-e994-48b7-aac7-d595ce81abd5
TableOfContents()

# ╔═╡ bfd62223-fdd7-47e9-8e53-0343ebe384db
md"""
## Load the data
"""

# ╔═╡ 88c7c602-9e9a-4def-8fd5-9ed4cfbafa6c
j, d = 3, 1

# ╔═╡ ed0c4e07-e0e7-4606-ab2d-1f4cd89ee20a
directory = "data/branch_j=$(j)_d=$(d)"

# ╔═╡ 69ecebd5-5e22-4fdc-b6ae-5282649a46eb
parameters, data_top, data_turn, data_bottom, data_connection_points =
    CGL.read_proof_witness(directory)

# ╔═╡ 72038085-72cd-4f45-8181-950632dc641b
md"""
## Check the data
"""

# ╔═╡ ddffe6ef-04e0-4d5d-b1dd-6f27c6029d48
CGL.check_proof_witness(
    parameters,
    data_top,
    data_turn,
    data_bottom,
    data_connection_points,
)

# ╔═╡ 0a4484f0-2b9a-4a12-893e-aa894c818c59
md"""
## Plot the data
"""

# ╔═╡ 73eff1b3-5cab-46bd-8371-b49c3d644a3f
br = CGL.CGLBranch.branch_epsilon(CGL.CGLBranch.sverak_initial(j, d)...)

# ╔═╡ 093f8c5a-11ba-49e0-b6f1-01d30736173e
md"""
### Full branch
"""

# ╔═╡ a3251cf9-03ca-4cb7-b371-c47dad1c3d00
pl_full = let
    pl = plot(xlabel = L"\epsilon", ylabel = L"\kappa", dpi = 300)

    plot!(
        pl,
        data_top.ϵ_lower,
        data_top.κ_exists,
        color = :green,
        linewidth = 5,
        label = "",
    )
    if !isnothing(data_turn)
        plot!(
            pl,
            data_turn.ϵ_exists,
            data_turn.κ_lower,
            color = :blue,
            linewidth = 5,
            label = "",
        )
    end
    if !isnothing(data_bottom)
        plot!(
            pl,
            data_bottom.ϵ_lower,
            data_bottom.κ_exists,
            color = :green,
            linewidth = 5,
            label = "",
        )
    end
    if !isempty(data_connection_points)
        scatter!(
            data_connection_points.ϵ,
            data_connection_points.κ_exists,
            markersize = 3,
            color = :red,
            label = "Connection points",
        )
    end

    # Fix limits
    #xlims!(pl, xlims(pl))
    #ylims!(pl, ylims(pl))

    plot!(pl, br.param, br.κ, color = :blue, linewidth = 1, label = "", linestyle = :dot)

    pl
end


# ╔═╡ ccdd845c-8c5f-4ca4-b4bc-a7e35954d645
md"""
### Beginning and end of branch
"""

# ╔═╡ 126de592-9cd0-493c-8ce8-70a9723989b2
pl_start = let data = data_top[1:100, :]
    ϵs = tuple.(data.ϵ_lower, data.ϵ_upper)
    κs_uniq = data.κ_uniq
    κs_exists = data.κ_exists

    pl = plot(xlabel = L"\epsilon", ylabel = L"\kappa", legend = :none)

    plot!(pl, vcat.(interval.(ϵs), interval.(κs_uniq)), linecolor = :red, color = :red)
    plot!(
        pl,
        vcat.(interval.(ϵs), interval.(κs_exists)),
        linecolor = :green,
        color = :green,
    )

    # Fix limits
    xlims!(pl, xlims(pl))
    ylims!(pl, ylims(pl))

    plot!(pl, br.param, br.κ, color = :blue, marker = :circle, markersize = 2)

    pl
end


# ╔═╡ ff0e261c-9ad3-4a42-96d3-ad45dbe8929a
pl_end = let

    data = if !isnothing(data_bottom)
        data_bottom[end-100:end, :]
    elseif !isnothing(data_turn)
        data_turn[end-100:end, :] # FIXME: Currently not handled
    else
        data_top[end-100:end, :]
    end
    ϵs = tuple.(data.ϵ_lower, data.ϵ_upper)
    κs_uniq = data.κ_uniq
    κs_exists = data.κ_exists

    pl = plot(xlabel = L"\epsilon", ylabel = L"\kappa", legend = :none)

    plot!(pl, vcat.(interval.(ϵs), interval.(κs_uniq)), linecolor = :red, color = :red)
    plot!(
        pl,
        vcat.(interval.(ϵs), interval.(κs_exists)),
        linecolor = :green,
        color = :green,
    )

    # Fix limits
    xlims!(pl, xlims(pl))
    ylims!(pl, ylims(pl))

    plot!(pl, br.param, br.κ, color = :blue, marker = :circle, markersize = 2)

    pl
end

# ╔═╡ 74711ce9-8317-4ee4-85d4-7317e5ea6d0c
md"""
### Connection points
"""

# ╔═╡ d7ab30ea-4813-4253-abc7-c27148962980
pl_top_turn = if isnothing(data_turn)
    missing
else
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

        ϵs_top = tuple.(data_top.ϵ_lower[end-1:end], data_top.ϵ_upper[end-1:end])
        κs_uniq_top = data_top.κ_uniq[end-1:end]
        κs_exists_top = data_top.κ_exists[end-1:end]

        ϵs_uniq_turn = data_turn.ϵ_uniq[i-1:i+1]
        ϵs_exists_turn = data_turn.ϵ_exists[i-1:i+1]
        κs_turn = tuple.(data_turn.κ_lower[i-1:i+1], data_turn.κ_upper[i-1:i+1])

        pl = plot(xlabel = L"\epsilon", ylabel = L"\kappa", legend = :none)

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
            data_connection_points.ϵ[1:1],
            data_connection_points.κ_exists[1:1],
            markersize = 3,
            color = :red,
        )

        xlims!(pl, xlims(pl))
        ylims!(pl, ylims(pl))

        plot!(pl, br.param, br.κ, color = :blue, marker = :circle, markersize = 2)
    end
end

# ╔═╡ e3f00871-77a2-4d2d-bf01-93076fdeca90
pl_bottom_turn = if isnothing(data_bottom)
    missing
else
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

        ϵs_bottom = tuple.(data_bottom.ϵ_lower[1:2], data_bottom.ϵ_upper[1:2])
        κs_uniq_bottom = data_bottom.κ_uniq[1:2]
        κs_exists_bottom = data_bottom.κ_exists[1:2]

        ϵs_uniq_turn = data_turn.ϵ_uniq[i-1:i+1]
        ϵs_exists_turn = data_turn.ϵ_exists[i-1:i+1]
        κs_turn = tuple.(data_turn.κ_lower[i-1:i+1], data_turn.κ_upper[i-1:i+1])

        pl = plot(xlabel = L"\epsilon", ylabel = L"\kappa", legend = :none)

        plot!(
            pl,
            vcat.(interval.(ϵs_exists_turn), interval.(κs_turn)),
            linecolor = :blue,
            color = :blue,
        )

        plot!(
            pl,
            vcat.(interval.(ϵs_bottom), interval.(κs_exists_bottom)),
            linecolor = :green,
            color = :green,
        )

        scatter!(
            data_connection_points.ϵ[2:2],
            data_connection_points.κ_exists[2:2],
            markersize = 3,
            color = :red,
        )

        xlims!(pl, xlims(pl))
        ylims!(pl, ylims(pl))

        plot!(pl, br.param, br.κ, color = :blue, marker = :circle, markersize = 2)
    end
end

# ╔═╡ 0c01333d-cd4d-4f46-89cc-e9de7fb94baa
md"""
### Other plots
"""

# ╔═╡ e14feb21-c51a-4add-af33-4a0684bae5f9
let
    categories = [
        "Existence top",
        "Existence turn",
        "Existence bottom",
        "Continuation top",
        "Continuation turn",
        "Continuation bottom",
    ]
    values = map(
        x -> isnan(x) ? 0.0 : x,
        [
            parameters.runtime_existence_top,
            parameters.runtime_existence_turn,
            parameters.runtime_existence_bottom,
            parameters.runtime_continuation_top,
            parameters.runtime_continuation_turn,
            parameters.runtime_continuation_bottom,
        ],
    )
    total_time = ceil(Int, sum(values))
    pie(categories, values, title = "Total computational time $(total_time ÷ 60) minutes")
end

# ╔═╡ 6b3c0193-2ff8-4f63-b6bb-be22f83ffdf1
let
    ϵ_width = if isnothing(data_turn)
        maximum(br.param)
    else
        maximum(data_turn.ϵ_exists)
    end
    κ_width = if isnothing(data_bottom)
        data_top.κ_exists[1] - br.κ[end]
    else
        data_top.κ_exists[1] - data_bottom.κ_exists[end]
    end

    pl = plot(yaxis = :log10)

    rad_top = data_top.ϵ_upper - data_top.ϵ_lower
    scatter!(
        pl,
        data_top.ϵ_lower,
        rad_top / ϵ_width,
        markerstrokewidth = 0,
        markersize = 2,
        label = "Top",
    )
    if !isnothing(data_turn)
        rad_turn = data_turn.κ_upper - data_turn.κ_lower
        scatter!(
            pl,
            data_turn.ϵ_exists,
            rad_turn / κ_width,
            markerstrokewidth = 0,
            markersize = 2,
            label = "Turn",
        )
    end
    if !isnothing(data_bottom)
        rad_bottom = data_bottom.ϵ_upper - data_bottom.ϵ_lower
        scatter!(
            pl,
            data_bottom.ϵ_lower,
            rad_bottom / ϵ_width,
            markerstrokewidth = 0,
            markersize = 2,
            label = "Bottom",
        )
    end
    pl
end

# ╔═╡ Cell order:
# ╠═2f550824-d714-11ee-1c46-9b93dd3db181
# ╠═ec89485a-4c4d-45d3-b2fe-1a1a1d61f147
# ╠═b073c9ef-0dcd-4496-b72b-f36fac091474
# ╠═ed29ab9a-e994-48b7-aac7-d595ce81abd5
# ╟─bfd62223-fdd7-47e9-8e53-0343ebe384db
# ╠═88c7c602-9e9a-4def-8fd5-9ed4cfbafa6c
# ╠═ed0c4e07-e0e7-4606-ab2d-1f4cd89ee20a
# ╠═69ecebd5-5e22-4fdc-b6ae-5282649a46eb
# ╟─72038085-72cd-4f45-8181-950632dc641b
# ╠═ddffe6ef-04e0-4d5d-b1dd-6f27c6029d48
# ╟─0a4484f0-2b9a-4a12-893e-aa894c818c59
# ╠═73eff1b3-5cab-46bd-8371-b49c3d644a3f
# ╟─093f8c5a-11ba-49e0-b6f1-01d30736173e
# ╟─a3251cf9-03ca-4cb7-b371-c47dad1c3d00
# ╟─ccdd845c-8c5f-4ca4-b4bc-a7e35954d645
# ╟─126de592-9cd0-493c-8ce8-70a9723989b2
# ╟─ff0e261c-9ad3-4a42-96d3-ad45dbe8929a
# ╟─74711ce9-8317-4ee4-85d4-7317e5ea6d0c
# ╟─d7ab30ea-4813-4253-abc7-c27148962980
# ╟─e3f00871-77a2-4d2d-bf01-93076fdeca90
# ╟─0c01333d-cd4d-4f46-89cc-e9de7fb94baa
# ╟─e14feb21-c51a-4add-af33-4a0684bae5f9
# ╟─6b3c0193-2ff8-4f63-b6bb-be22f83ffdf1
