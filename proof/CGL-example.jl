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

# ╔═╡ 00d28bd0-7b41-11ef-045e-5d50ae3a6149
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using Arblib
    using CGL
    using DataFrames
    using IntervalArithmetic: interval
    using LaTeXStrings
    using Plots
    using PlutoUI
    using Printf
    using StaticArrays

    # Set the bits of precision used for the computations with Arblib
    setprecision(Arb, 128)

    ENV["LD_LIBRARY_PATH"] = "/home/joeldahne/Programs/capd/lib/"

    nothing
end

# ╔═╡ e16c81ab-cc21-4008-b170-b36d6439fea5
md"""
# Detail example for CGL
In this notebook we go through a detailed example with the computations that are required to prove the existence of a branch of self-similar singular solutions to the complex Ginzburg-Landau equation, as well as to determine the number of critical points of the profiles along the branch.

The notebook corresponds to Section 6.1 in the paper with which this repository is associated, it follows the same structure, and all numbers in that section of the paper are coming from this notebook. The notebook serves to makes it explicit what code is run to get these results. Since the output of the notebook is used in the paper, part of the notebook is concerned with formatting the output in a way that is suitable for inclusion in LaTeX.
"""

# ╔═╡ 07046b42-2935-4eb3-937e-412a6e58e309
TableOfContents()

# ╔═╡ 5c38a867-e6fd-4223-85f4-3c8445e7cdff
md"""
We go through the details for the first branch in *Case I*, corresponding to $d = 1$, $\sigma = 2.3$ and $\delta = 0$. In the code this is given by $j = 1$ and $d = 1$.
"""

# ╔═╡ 171e328b-7d09-447c-b211-3fe17dbfdb92
j, d = 1, 1

# ╔═╡ 445f4f90-b01b-4772-ac8d-f643d14fbf33
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save CheckBox(default = false))
"""

# ╔═╡ cb45ecca-be99-403e-ad5e-8eefbb4bd2d0
md"""
## Numerical approximation of the branch
The first step is to compute a numerical approximation of the branch. The heavy lifting for handling this is done by the package [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl), which does automatic continuation of solutions. The implementation for CGL is handled by the function `CGL.CGLBranch.branch_epsilon`, which returns a representation of the branch.
"""

# ╔═╡ 69dce4dd-7ba2-4264-b639-307b66924820
branch_approximation = CGL.CGLBranch.branch_epsilon(CGL.CGLBranch.sverak_initial(j, d)...)

# ╔═╡ 6fb987e8-8731-4886-9af6-4f6958ed03a3
md"""
Each point on the curve is represented by a tripple $(\epsilon, \mu, \kappa)$. Note that the parameter $\gamma$ (used to parameterize the manifold at infinity) does not directly appear, it is only used as an internal variable at this stage. We can plot the $(\epsilon, \kappa)$-projection of this curve,
"""

# ╔═╡ 5638715e-2f69-45ba-b45b-5739087d1b04
let pl = plot(legend = :none, xlabel = L"\epsilon", ylabel = L"\kappa")
    plot!(pl, branch_approximation.param, branch_approximation.κ, linestyle = :dash)

    save && savefig(pl, "figures/CGL-example-approximation-kappa.pdf")

    pl
end

# ╔═╡ 5646414a-75b2-4e11-a313-e09a46e6d603
md"""
as well as the $(\epsilon, \mu)$ projection
"""

# ╔═╡ 7b19b234-d906-425a-80f1-cd7f35c93a3b
let pl = plot(legend = :none, xlabel = L"\epsilon", ylabel = L"\mu")
    plot!(pl, branch_approximation.param, branch_approximation.µ, linestyle = :dash)

    save && savefig(pl, "figures/CGL-example-approximation-mu.pdf")

    pl
end

# ╔═╡ d27b3ed4-aeb6-486a-8887-3918df1b7ee7
md"""
## Verification of the beginning of the branch
Verifying the full branch is very computationally expensive, for that reason we start by verifying only a small part at the beginning of the branch. The verification of the full branch is discussed in the end of this notebook, and is based on precomputed data. The first two points on the branch in our numerical approximation are given by
"""

# ╔═╡ 1694a993-09ca-49b4-8de1-09c7d3f14ba6
point₁ = (
    ϵ = branch_approximation.param[1],
    µ = branch_approximation.μ[1],
    κ = branch_approximation.κ[1],
)

# ╔═╡ cb8fb346-de8f-4a3d-9cb6-d9553054391c
point₂ = (
    ϵ = branch_approximation.param[2],
    µ = branch_approximation.μ[2],
    κ = branch_approximation.κ[2],
)

# ╔═╡ 6b5aaff9-77b0-4e50-84f0-ad97ea43b875
ϵ₁ = point₁.ϵ

# ╔═╡ 8693456d-274e-499f-8311-cc95dac1f46a
ϵ₂ = point₂.ϵ

# ╔═╡ 468bdbe1-b0f5-4742-94d5-24cd9cefd798
md"""
Our goal will be to prove the existence of a continuous curve of solutions, starting at $\epsilon_1$, and ending at $\epsilon_2$. This proof will come together as a collection of boxes covering the interval, where each box is proved to contain a solution for each $\epsilon$ in the box. The picture we will end up with will look like this:
"""

# ╔═╡ bac3231c-e121-4468-8c27-bd1837ba8648
md"""
Recall that the goal is to prove that there exists a continuous curve, $(\epsilon(s), \mu(s), \kappa(s))$, such that for each point on this curve, the CGL equation has a self-similar singular solution. The $(\epsilon, \kappa)$-projection of this curve will be contained inside the boxes in the above figure. For the beginning of the branch, which is what we are considering now, we can take the curve to be parameterized by $\epsilon$. We are then looking for a curve of the form $(\epsilon, \mu(\epsilon), \kappa(\epsilon))$, with $\epsilon \in [\epsilon_1, \epsilon_2]$. In practice, the curve we compute will also contain the value of $\gamma$, so the curve is given by $(\epsilon, \mu(\epsilon), \gamma(\epsilon), \kappa(\epsilon))$, But in the statement of the theorem, the $\gamma$ is dropped.
"""

# ╔═╡ 1a0ee2d2-2bd0-49b3-85c6-ab37424dc562
md"""
Before verifying that the curve segments join together, let us however first take a look at the computation of the boxes $\boldsymbol{\mu}_{i}^{(e)} \times \boldsymbol{\gamma}_{i}^{(e)} \times \boldsymbol{\kappa}_{i}^{(e)}$ and $\boldsymbol{\mu}_{i}^{(u)} \times \boldsymbol{\gamma}_{i}^{(u)} \times \boldsymbol{\kappa}_{i}^{(u)}$. This computation is handled by the function `branch_segment_existence_fix_epsilon`. It takes as argument the $\mu$, $\kappa$ and $\epsilon$ coefficients for two points on the curve, and then tries to split the part between these two points into boxes where existence and uniqueness can be verified.

The splitting is based on a dynamic bisection approach. It starts with the full interval $[\epsilon_1, \epsilon_2]$, and if that fails it bisects this interval into two smaller ones, the procedure is then recursively iterated for two smaller intervals.
"""

# ╔═╡ c9565b23-bf18-4f71-907a-b53bdaa00a53
ξ₁, λ = CGL.sverak_params(Arb, j, d, ξ₁_for_branch = true)[end-1:end]

# ╔═╡ d6a4a3bb-fdae-4e2b-a2c0-47ab4ab805f7
ϵs, exists, uniqs = CGL.branch_segment_existence_fix_epsilon(
    (Arb(point₁.μ), Arb(point₂.μ)),
    (Arb(point₁.κ), Arb(point₂.κ)),
    (Arf(point₁.ϵ), Arf(point₂.ϵ)),
    ξ₁,
    λ,
    depth_start = 0,
    verbose = true,
)[1:3]

# ╔═╡ 1d47924a-05c1-4d43-8c2d-649159e9c714
let pl = plot(xlabel = L"\epsilon", ylabel = L"\kappa")
    plot!(
        pl,
        vcat.(interval.(ϵs), interval.(getindex.(exists, 4))),
        label = "Enclosure of branch",
    )

    scatter!(pl, [point₁.ϵ, point₂.ϵ], [point₁.κ, point₂.κ], label = "Approximations")

    pl
end

# ╔═╡ 97dc6eba-2235-48a3-9829-f65ee4f8833a
N = length(ϵs)

# ╔═╡ fdffe508-4a94-4e92-b524-a58db5049d2a
md"""
The boxes split the interval $[\epsilon_1, \epsilon_2]$ into $N =$ $N subintervals,

$$[\epsilon_1, \epsilon_2] = \bigcup_{i = 1}^{N} \boldsymbol{\epsilon}_{i}.$$

For each subinterval, $\boldsymbol{\epsilon}_{i}$, the proof of existence, and uniqueness, follows exactly the same approach as in the example for NLS, see the `NLS_example.jl` notebook. The only difference is that for NLS we fixed $\epsilon = 0$, whereas in this case we do the computations with $\epsilon$ represented by the interval $\boldsymbol{\epsilon}_i$. The properties of interval arithmetic then ensures us that any results we get out of this procedure, in this case enclosures of existence and uniqueness, are valid for any $\epsilon \in \boldsymbol{\epsilon}_i$.

More precisely, the procedure produces two boxes,

$$\boldsymbol{\mu}_{i}^{(e)} \times \boldsymbol{\gamma}_{i}^{(e)} \times \boldsymbol{\kappa}_{i}^{(e)} \subsetneq
\boldsymbol{\mu}_{i}^{(u)} \times \boldsymbol{\gamma}_{i}^{(u)} \times \boldsymbol{\kappa}_{i}^{(u)},$$

where for all $\epsilon \in \boldsymbol{\epsilon}_{i}$ the function $G$
is proved to have a zero in the inner box, and this zero is unique in
the outer box. Moreover, the Jacobian of $G$ is guaranteed to be
invertible for all $\epsilon \in \boldsymbol{\epsilon}_{i}$, and
parameters in the outer box. Since $G$ is continuously
differentiable with respect to all of its arguments, it follows by the
implicit function theorem that there exists a continuous curve of
solutions

$$(\epsilon, \mu_{i}(\epsilon), \gamma_{i}(\epsilon), \kappa_{i}(\epsilon)),$$

defined for $\epsilon \in \boldsymbol{\epsilon}_{i}$.

This gives us $N$ segments of curves, to prove that we have a continuous curve on the whole interval $[\epsilon_1, \epsilon_2]$ we need to verify that these $N$ segments in fact join together to form one long curve. 
"""

# ╔═╡ 7ebc1c49-969b-44e2-a236-294f3c4e3231
all(exist -> all(isfinite, exist), exists) # Check that is succeeded

# ╔═╡ 46df8e33-82ba-41c4-842e-b167f03e0220
md"""
For the return value `ϵs`, the value at index `i`, i.e. `ϵs[i]`, corresponds to the interval $\boldsymbol{\epsilon}_i$. The `exists` and `uniqs` results gives the boxes for existence and uniquenes respectively. The value of $\boldsymbol{\epsilon}_i$ corresponds to
"""

# ╔═╡ 76a026ad-60fc-42c7-91c1-ca51cf8bcbde
ϵs[1]

# ╔═╡ d8bfb93a-8b10-4ffa-8d8f-afbee26407b0
md"""
and the box $\boldsymbol{\mu}_{1}^{(e)} \times \boldsymbol{\gamma}_{1}^{(e)} \times \boldsymbol{\kappa}_{1}^{(e)}$ would correspond to
"""

# ╔═╡ 4efd50df-6b92-475b-abaa-740ecc39fe57
exists[1]

# ╔═╡ 39308e62-fb03-4571-ac79-4a0e0f87bc90
md"""
where `exists[1][1]` is $\boldsymbol{\mu}_{1}^{(e)}$, `exists[1][2]` and `exists[1][3]` corresponds to the real and imaginary parts of $\boldsymbol{\gamma}_{1}^{(e)}$ and `exists[1][4]` gives $\boldsymbol{\kappa}_{1}^{(e)}$. Similarly, the box $\boldsymbol{\mu}_{1}^{(u)} \times \boldsymbol{\gamma}_{1}^{(u)} \times \boldsymbol{\kappa}_{1}^{(u)}$ corresponds to
"""

# ╔═╡ 3c542720-542a-4ff7-8ca4-c620582d95ab
uniqs[1]

# ╔═╡ 17cb9194-be0d-4596-b761-4d5747f759d3
md"""
Plotting the $(\epsilon, \kappa)$-projection of the boxes for existence and uniqueness gives us:
"""

# ╔═╡ 0fadf749-fc85-4223-925f-bf5b4b88239e
let pl = plot(xlabel = L"\epsilon", ylabel = L"\kappa")
    plot!(pl, vcat.(interval.(ϵs), interval.(getindex.(uniqs, 4))), label = "Uniqueness")

    plot!(pl, vcat.(interval.(ϵs), interval.(getindex.(exists, 4))), label = "Existence")

    scatter!(pl, [point₁.ϵ, point₂.ϵ], [point₁.κ, point₂.κ], label = "Approximations")

    save && savefig(pl, "figures/CGL-example-beginning.pdf")

    pl
end

# ╔═╡ 322ab6ec-58c2-46f7-85d9-7d59db7b8fd3
md"""
Now, let us come back to the question of whether the $N$ separate curve segments can be joined together into one continuous curve. What one needs to avoid is the situation depicted in the figure below, where we have two separate curves and our enclosures go from enclosing the bottom curve, to enclosing the upper curve. To ensure that this situation does not happen it is enough to verify that

$$\boldsymbol{\mu}_{i}^{(e)} \times \boldsymbol{\gamma}_{i}^{(e)} \times \boldsymbol{\kappa}_{i}^{(e)}
  \subseteq
  \boldsymbol{\mu}_{i - 1}^{(u)} \times \boldsymbol{\gamma}_{i - 1}^{(u)} \times \boldsymbol{\kappa}_{i - 1}^{(u)}$$

holds for all $i = 2, \dots, N$, i.e., the enclosure of existence for the current interval is included in the enclosure of uniqueness for the previous interval. Note that in the figure below this condition is not satisfied for the middle parts.
"""

# ╔═╡ 3749abf9-e8aa-4761-8d10-004f0e50976e
let pl = plot(axis = ([], false), grid = false)
    xs = range(-0.1, 0.1, 1000)
    intervals = CGL.mince(Arb((-0.1, 0.1)), 16)

    f_1 = x -> 0.998cos(x)
    ys_1 = f_1.(xs)
    ys_exists_1 = add_error.(f_1.(intervals) .+ 0.0004, Mag(0.0008))
    ys_uniqs_1 = add_error.(f_1.(intervals) .+ 0.0004, Mag(0.0014))

    f_2 = x -> 2 - cos(x)
    ys_2 = f_2.(xs)
    ys_exists_2 = add_error.(f_2.(intervals) .- 0.0004, Mag(0.0008))
    ys_uniqs_2 = add_error.(f_2.(intervals) .- 0.0004, Mag(0.0014))

    plot!(pl, xs, ys_1, color = :blue, linewidth = 3, label = "")
    plot!(pl, xs, ys_2, color = :blue, linewidth = 3, label = "")

    plot!(
        pl,
        vcat.(interval.(intervals), interval.(ys_uniqs_1))[1:8],
        color = :yellow,
        label = "Uniqueness",
    )
    plot!(
        pl,
        vcat.(interval.(intervals), interval.(ys_exists_1))[1:8],
        color = :green,
        label = "Existence",
    )

    plot!(
        pl,
        vcat.(interval.(intervals), interval.(ys_uniqs_2))[9:16],
        color = :yellow,
        label = "",
    )
    plot!(
        pl,
        vcat.(interval.(intervals), interval.(ys_exists_2))[9:16],
        color = :green,
        label = "",
    )

    save && savefig(pl, "figures/CGL-example-uniqueness-issue.pdf")

    pl
end

# ╔═╡ f8822ca2-78b1-4ef0-ab01-ac4d0dac20f7
md"""
In our case, this condition is however easily verified to be true for all parts.
"""

# ╔═╡ 43d2ef90-83cf-4497-93ca-37904f8ca3d6
all(eachindex(exists, uniqs)[2:end]) do i
    all(Arblib.contains.(uniqs[i-1], exists[i]))
end

# ╔═╡ e669a71f-393a-4466-893c-f923b240c34f
md"""
In some cases it will however not be possible to directly verify this condition, this happens further down on some of the branches. In that case we divide the intervals $\boldsymbol{\epsilon}_i$ into smaller subintervals, this gives tighter enclosures for the existence, which makes the condition easier to verify.
"""

# ╔═╡ bdf81d6d-9a5d-46e6-83d1-4169265eef65
md"""
The last remaining part is to count the number of critical points of the profile along the curve. This is done for each box independently, and the procedure is the same as in the example for NLS, see the NLS-example.jl notebook.
"""

# ╔═╡ 768c3409-40ca-4e69-bd1c-eb622bbffa5a
num_critical_points = map(ϵs, exists) do ϵ, exist
    μ, γ, κ = exist[1], Acb(exist[2], exist[3]), exist[4]
    success, critical_points, verified =
        CGL.count_critical_points(μ, γ, κ, Arb(ϵ), ξ₁, λ)
    success ? length(critical_points) : missing
end

# ╔═╡ 8518d0f0-6790-4cab-8611-801afd50413b
md"""
We can verify that indeed all boxes reported having exactly $j - 1$ critical points in this case.
"""

# ╔═╡ e4d3f0ea-e4bd-475e-8edc-829bdbbefba7
all(==(j - 1), num_critical_points)

# ╔═╡ f8d78162-53c2-4750-81c6-4db5b8a67f28
md"""
In this case the process succeeded for all subintervals, but in some cases further down the branch the error bounds are larger and the process sometimes fails. In that case we split each box into smaller pieces, allowing us to get tigher enclosures, and increasing the chance of success. This extra splitting comes with a computational cost, but does allow us to count the number of critical points along the whole branch.
"""

# ╔═╡ bbde7e0b-517d-4154-89f4-df5d1a409a94
md"""
## Verification of full branch
The verification of the full branch follows the same procedure as above, there are however a few more things that needs to be taken into account.

At the beginning of the branch it was possible to parameterize the curve using $\epsilon$, once we reach the turning point of the branch this is however no longer possible, and a different parameterization has to be used. For this we split the curve into three different parts, one top part where the parameterization is done in $\epsilon$, one turn part where the parameterization is done in $\kappa$ and one bottom part where the parameterization is again done in $\epsilon$. The top, turn and bottom parts of the branch are shown in the figure below, top and bottom in green and the turn in blue.
"""

# ╔═╡ 951fe56e-9883-49f8-8845-34f137a51bd4
md"""
For the top and bottom parts the procedure is exactly as described above. For the turning part we have to flip the role of $\epsilon$ and $\kappa$, instead of splitting it into intervals in $\epsilon$, we split it in $\kappa$ as

$$\bigcup_{i = 1}^{N} \boldsymbol{\kappa}_{i}.$$

And when appying the interval Newton method on $G$, we now fix $\kappa$ and look for a zero in $\mu$, $\gamma$ and $\epsilon$ instead. The main difference this gives, is that we have to compute the Jacobian of $G$ with respect to $\mu$, $\gamma$ and $\epsilon$, instead of with respect to $\mu$, $\gamma$ and $\kappa$. Other than that, the procedure is the same.

With this we get three continuous curves, one for the top part, one for the turning part and one for the bottom part. The only thing that remains is to prove that these three curves can be joined together into one curve. In this case we can't use the same approach as described above, that requires the curves to be parameterized by the same variable. Instead we will prove that the curves can be connected by finding a point that is proved to be contained in both curves, in which the two curves has to be joined together.

This process, for joining the top part and the turning part, is depicted in the figure below. The green boxes represent the uniqueness for the top part and the blue boxes the uniqueness for the turning part, the red point is very accurately computed enclosure of a solution. Since the red point lies within the box of uniqueness for the top part, it has to be part of the top curve. Since it also lies within the box of uniqueness for the turning part, it has to also be part of the turning curve. That means that these two curves share a common points, and must hence be the same.
"""

# ╔═╡ b2be1faf-d2a0-4a5e-9546-ffbadcdd753a
md"""
For connecting the turning part and the bottom part the procedure is the same, and depicted in the figure below.
"""

# ╔═╡ 7608c6e2-ecb7-4c19-b60b-a4235ee017dc
md"""
## Proof data
"""

# ╔═╡ 95ae3d57-3e3d-4d37-a0c5-7feeb655341f
md"""
The data on which the proof for the full branch is based on can be loaded using the function `read_proof_witness`.
"""

# ╔═╡ 72f41478-3d5a-47cd-a4a3-e9aa8d0f1672
parameters, data_top, data_turn, data_bottom, data_connection_points =
    CGL.read_proof_witness("data/branch_j=$(j)_d=$(d)/");

# ╔═╡ a7e3d3d2-7cd0-4fe3-afe4-dd02499d27c9
let
    pl = plot(legend = :none, xlabel = L"\epsilon", ylabel = L"\kappa")

    top_indices = round.(Int, range(1, nrow(data_top), length = 1000))
    turn_indices = round.(Int, range(1, nrow(data_turn), length = 1000))
    bottom_indices = round.(Int, range(1, nrow(data_bottom), length = 1000))

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

    save && savefig(pl, "figures/CGL-example-parts.pdf")

    pl
end

# ╔═╡ 5d50bb60-17ab-43ad-a20c-f16cd595dc20
let pl = plot(xlabel = L"\epsilon", ylabel = L"\kappa", legend = :none)
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

    boxes_top =
        vcat.(interval.(data_top.ϵ_lower, data_top.ϵ_upper), interval.(data_top.κ_uniq))

    boxes_turn =
        vcat.(interval.(data_turn.ϵ_uniq), interval.(data_turn.κ_lower, data_turn.κ_upper))

    # Convenient way of picking good xlims and ylims values
    xlims_value, ylims_value = let pl_tmp = plot()
        plot!(pl_tmp, boxes_turn[i-1:i+1])

        plot!(pl_tmp, boxes_top[end-1:end])

        xlims(pl_tmp), ylims(pl_tmp)
    end

    xlims!(pl, xlims_value)
    ylims!(pl, ylims_value)

    plot!(pl, boxes_turn, linecolor = :blue, color = :blue)

    plot!(pl, boxes_top, linecolor = :green, color = :green)

    scatter!(
        data_connection_points.ϵ[1:1],
        data_connection_points.κ_exists[1:1],
        markersize = 3,
        color = :red,
    )

    save && savefig(pl, "figures/CGL-example-top-turn.pdf")

    pl
end

# ╔═╡ 18b4e466-3631-4014-a6bb-ee9a69ebb337
let pl = plot(xlabel = L"\epsilon", ylabel = L"\kappa", legend = :none)
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

    boxes_bottom =
        vcat.(
            interval.(data_bottom.ϵ_lower, data_bottom.ϵ_upper),
            interval.(data_bottom.κ_uniq),
        )

    boxes_turn =
        vcat.(interval.(data_turn.ϵ_uniq), interval.(data_turn.κ_lower, data_turn.κ_upper))

    # Convenient way of picking good xlims and ylims values
    xlims_value, ylims_value = let pl_tmp = plot()
        plot!(pl_tmp, boxes_turn[i-1:i+1])

        plot!(pl_tmp, boxes_bottom[1:2])

        xlims(pl_tmp), ylims(pl_tmp)
    end

    xlims!(pl, xlims_value)
    ylims!(pl, ylims_value)

    plot!(pl, boxes_turn, linecolor = :blue, color = :blue)

    plot!(pl, boxes_bottom, linecolor = :green, color = :green)

    scatter!(
        data_connection_points.ϵ[2:2],
        data_connection_points.κ_exists[2:2],
        markersize = 3,
        color = :red,
    )

    save && savefig(pl, "figures/CGL-example-turn-bottom.pdf")

    pl
end

# ╔═╡ 9654fafc-fcd4-4f70-b552-7e3ef25142a0
md"""
The function `check_proof_witness` can be used to do a sanity check of the data. Among a few other things, it checks that the subintervals follow consecutively after each other, that the computed enclosures are finite and that the existence and uniqueness enclosures satisfy the required inclusion for continuation. All of these properties are guaranteed by the way that the results are computed, but an extra check reduces the risk of any mistakes.
"""

# ╔═╡ 87494e0d-c7b5-4027-8d1b-a66e48f19750
CGL.check_proof_witness(
    parameters,
    data_top,
    data_turn,
    data_bottom,
    data_connection_points,
)

# ╔═╡ 2d651446-c65c-46ba-a51f-4fdcb05a60f3
md"""
We can also get from statistics from the data, such as the number of subintervals for the top, turn and bottom parts.
"""

# ╔═╡ b49b3fb8-d281-49b3-842e-cad10a3d1e7a
nrow(data_top), nrow(data_turn), nrow(data_bottom)

# ╔═╡ c9bfb637-7dc5-4c0b-82b3-4e8f38b534c8
md"""
The runtime (in seconds) for proving the existence of the different parts:
"""

# ╔═╡ df33530f-1cbf-4702-b10b-03c2a55558cb
runtime_top =
    round(Int, parameters.runtime_existence_top + parameters.runtime_continuation_top)

# ╔═╡ 4ce9d880-1844-431e-be38-859450b721ac
runtime_turn = parameters.runtime_existence_turn + parameters.runtime_continuation_turn

# ╔═╡ 2cb2f265-7029-4df5-aca2-2261d890e6eb
runtime_bottom =
    parameters.runtime_existence_bottom + parameters.runtime_continuation_bottom

# ╔═╡ 182253a0-771d-4b60-8241-8bcaf1dc9b77
md"""
The runtime (in seconds) for counting the number of critical points for the different parts.
"""

# ╔═╡ 0ad117ee-c0d2-4bdb-9dc0-ebec8c46a394
runtime_critical_points_top = runtime_top

# ╔═╡ 454407dc-ba50-475c-b555-cff46f21b405
runtime_critical_points_turn = runtime_turn

# ╔═╡ 168cdbd8-7125-4dbb-befc-06b5dd05e888
runtime_critical_points_bottom = runtime_bottom

# ╔═╡ ea173928-4f1c-44b2-b62a-bdc67380c08e
md"""
We can get the corresponding number of core hours by scaling by the number of cores of the system used
"""

# ╔═╡ 2a4720ed-d503-4e09-9276-ed96f01b27f1
num_cores = 256

# ╔═╡ 91fb624a-a91e-4646-8eb9-1e298220d153
core_hours_top = runtime_top * num_cores / 3600

# ╔═╡ f13e187c-14f9-4302-b621-17b6f500f7d4
core_hours_turn = runtime_turn * num_cores / 3600

# ╔═╡ e07ba4fa-8e6f-41d1-bbe9-2c0924519e6d
core_hours_bottom = runtime_bottom * num_cores / 3600

# ╔═╡ f418cbcc-c6c5-414d-a364-04d8139c72fd
core_hours_critical_points_top = runtime_critical_points_top * num_cores / 3600

# ╔═╡ 65e32152-f61b-4b45-9483-248d6dd894e0
core_hours_critical_points_turn = runtime_critical_points_turn * num_cores / 3600

# ╔═╡ ea2169a6-3089-4c29-b8e5-0029acae31c6
core_hours_critical_points_bottom = runtime_critical_points_bottom * num_cores / 3600

# ╔═╡ dd21c22a-ae05-4878-8d3d-a52526bd9dcd
md"""
## Formatting output for paper
"""

# ╔═╡ 69d5fb8b-1d74-4f29-8621-8d5b06172e4e
point₁_latex = """
\\begin{equation*}
  \\epsilon_1 = $(point₁.ϵ),\\quad
  \\mu_1 = $(point₁.μ),\\quad
  \\kappa_1 = $(point₁.κ)
\\end{equation*}
"""

# ╔═╡ 32ffd546-b800-4eed-972a-a579b220b01f
point₂_latex = """
\\begin{equation*}
  \\epsilon_2 = $(point₂.ϵ),\\quad
  \\mu_2 = $(point₂.μ),\\quad
  \\kappa_2 = $(point₂.κ).
\\end{equation*}
"""

# ╔═╡ 5556da81-377e-4130-ba12-291a2e5c231a
md"""
The approximations for the first two points on our curve are given by
"""

# ╔═╡ 3a9e1303-b709-4e5a-a9d3-8bdf9af24c0e
latexstring(point₁_latex)

# ╔═╡ a7860708-b51c-4a9f-ac30-9f71b683528a
md"""
and
"""

# ╔═╡ d6ac9d6f-c589-4264-94d2-5b519b83fa80
latexstring(point₂_latex)

# ╔═╡ d5b23ed8-eef8-4d00-9e8a-73caa43a4737
print(point₁_latex)

# ╔═╡ cd4429fa-28d6-4c02-b585-8cd49af2d281
print(point₂_latex)

# ╔═╡ c2e72d74-c75d-40a4-ac5f-458c143f61cd
table_latex = """
Subintervals & $(nrow(data_top)) & $(nrow(data_turn)) & $(nrow(data_bottom))\\\\
Runtime existence of curve (seconds) & $(round(Int, runtime_top)) & $(round(Int, runtime_turn)) & $(round(Int, runtime_bottom))\\\\
Runtime existence of curve (core hours) & $(round(Int, core_hours_top)) & $(round(Int, core_hours_turn)) & $(round(Int, core_hours_bottom))\\\\
Runtime critical points (seconds) & $(round(Int, runtime_critical_points_top)) & $(round(Int, runtime_critical_points_turn)) & $(round(Int, runtime_critical_points_bottom))\\\\
Runtime critical points (core hours) & $(round(Int, core_hours_critical_points_top)) & $(round(Int, core_hours_critical_points_turn)) & $(round(Int, core_hours_critical_points_bottom))\\\\
"""

# ╔═╡ 8e79059c-bea8-47ae-94d2-78b6d2f4b390
print(table_latex)

# ╔═╡ Cell order:
# ╟─e16c81ab-cc21-4008-b170-b36d6439fea5
# ╠═00d28bd0-7b41-11ef-045e-5d50ae3a6149
# ╠═07046b42-2935-4eb3-937e-412a6e58e309
# ╟─5c38a867-e6fd-4223-85f4-3c8445e7cdff
# ╠═171e328b-7d09-447c-b211-3fe17dbfdb92
# ╟─445f4f90-b01b-4772-ac8d-f643d14fbf33
# ╟─cb45ecca-be99-403e-ad5e-8eefbb4bd2d0
# ╠═69dce4dd-7ba2-4264-b639-307b66924820
# ╟─6fb987e8-8731-4886-9af6-4f6958ed03a3
# ╟─5638715e-2f69-45ba-b45b-5739087d1b04
# ╟─5646414a-75b2-4e11-a313-e09a46e6d603
# ╟─7b19b234-d906-425a-80f1-cd7f35c93a3b
# ╟─d27b3ed4-aeb6-486a-8887-3918df1b7ee7
# ╠═1694a993-09ca-49b4-8de1-09c7d3f14ba6
# ╠═cb8fb346-de8f-4a3d-9cb6-d9553054391c
# ╠═6b5aaff9-77b0-4e50-84f0-ad97ea43b875
# ╠═8693456d-274e-499f-8311-cc95dac1f46a
# ╟─468bdbe1-b0f5-4742-94d5-24cd9cefd798
# ╟─1d47924a-05c1-4d43-8c2d-649159e9c714
# ╟─bac3231c-e121-4468-8c27-bd1837ba8648
# ╠═97dc6eba-2235-48a3-9829-f65ee4f8833a
# ╟─fdffe508-4a94-4e92-b524-a58db5049d2a
# ╟─1a0ee2d2-2bd0-49b3-85c6-ab37424dc562
# ╠═c9565b23-bf18-4f71-907a-b53bdaa00a53
# ╠═d6a4a3bb-fdae-4e2b-a2c0-47ab4ab805f7
# ╠═7ebc1c49-969b-44e2-a236-294f3c4e3231
# ╟─46df8e33-82ba-41c4-842e-b167f03e0220
# ╠═76a026ad-60fc-42c7-91c1-ca51cf8bcbde
# ╟─d8bfb93a-8b10-4ffa-8d8f-afbee26407b0
# ╠═4efd50df-6b92-475b-abaa-740ecc39fe57
# ╟─39308e62-fb03-4571-ac79-4a0e0f87bc90
# ╠═3c542720-542a-4ff7-8ca4-c620582d95ab
# ╟─17cb9194-be0d-4596-b761-4d5747f759d3
# ╟─0fadf749-fc85-4223-925f-bf5b4b88239e
# ╟─322ab6ec-58c2-46f7-85d9-7d59db7b8fd3
# ╟─3749abf9-e8aa-4761-8d10-004f0e50976e
# ╟─f8822ca2-78b1-4ef0-ab01-ac4d0dac20f7
# ╠═43d2ef90-83cf-4497-93ca-37904f8ca3d6
# ╟─e669a71f-393a-4466-893c-f923b240c34f
# ╟─bdf81d6d-9a5d-46e6-83d1-4169265eef65
# ╠═768c3409-40ca-4e69-bd1c-eb622bbffa5a
# ╟─8518d0f0-6790-4cab-8611-801afd50413b
# ╠═e4d3f0ea-e4bd-475e-8edc-829bdbbefba7
# ╟─f8d78162-53c2-4750-81c6-4db5b8a67f28
# ╟─bbde7e0b-517d-4154-89f4-df5d1a409a94
# ╟─a7e3d3d2-7cd0-4fe3-afe4-dd02499d27c9
# ╟─951fe56e-9883-49f8-8845-34f137a51bd4
# ╟─5d50bb60-17ab-43ad-a20c-f16cd595dc20
# ╟─b2be1faf-d2a0-4a5e-9546-ffbadcdd753a
# ╟─18b4e466-3631-4014-a6bb-ee9a69ebb337
# ╟─7608c6e2-ecb7-4c19-b60b-a4235ee017dc
# ╟─95ae3d57-3e3d-4d37-a0c5-7feeb655341f
# ╠═72f41478-3d5a-47cd-a4a3-e9aa8d0f1672
# ╟─9654fafc-fcd4-4f70-b552-7e3ef25142a0
# ╠═87494e0d-c7b5-4027-8d1b-a66e48f19750
# ╟─2d651446-c65c-46ba-a51f-4fdcb05a60f3
# ╠═b49b3fb8-d281-49b3-842e-cad10a3d1e7a
# ╟─c9bfb637-7dc5-4c0b-82b3-4e8f38b534c8
# ╠═df33530f-1cbf-4702-b10b-03c2a55558cb
# ╠═4ce9d880-1844-431e-be38-859450b721ac
# ╠═2cb2f265-7029-4df5-aca2-2261d890e6eb
# ╟─182253a0-771d-4b60-8241-8bcaf1dc9b77
# ╠═0ad117ee-c0d2-4bdb-9dc0-ebec8c46a394
# ╠═454407dc-ba50-475c-b555-cff46f21b405
# ╠═168cdbd8-7125-4dbb-befc-06b5dd05e888
# ╟─ea173928-4f1c-44b2-b62a-bdc67380c08e
# ╠═2a4720ed-d503-4e09-9276-ed96f01b27f1
# ╠═91fb624a-a91e-4646-8eb9-1e298220d153
# ╠═f13e187c-14f9-4302-b621-17b6f500f7d4
# ╠═e07ba4fa-8e6f-41d1-bbe9-2c0924519e6d
# ╠═f418cbcc-c6c5-414d-a364-04d8139c72fd
# ╠═65e32152-f61b-4b45-9483-248d6dd894e0
# ╠═ea2169a6-3089-4c29-b8e5-0029acae31c6
# ╟─dd21c22a-ae05-4878-8d3d-a52526bd9dcd
# ╟─69d5fb8b-1d74-4f29-8621-8d5b06172e4e
# ╟─32ffd546-b800-4eed-972a-a579b220b01f
# ╟─5556da81-377e-4130-ba12-291a2e5c231a
# ╟─3a9e1303-b709-4e5a-a9d3-8bdf9af24c0e
# ╟─a7860708-b51c-4a9f-ac30-9f71b683528a
# ╟─d6ac9d6f-c589-4264-94d2-5b519b83fa80
# ╠═d5b23ed8-eef8-4d00-9e8a-73caa43a4737
# ╠═cd4429fa-28d6-4c02-b585-8cd49af2d281
# ╟─c2e72d74-c75d-40a4-ac5f-458c143f61cd
# ╠═8e79059c-bea8-47ae-94d2-78b6d2f4b390
