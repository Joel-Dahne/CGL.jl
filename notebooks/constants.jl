### A Pluto.jl notebook ###
# v0.19.22

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

# ╔═╡ 8d7322fe-e018-11ed-1579-7daddd6ce28f
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using Arblib
    using DifferentialEquations
    using Folds
    using GinzburgLandauSelfSimilarSingular
    using LaTeXStrings
    using NLsolve
    using Plots
    using Plots.PlotMeasures
    using PlutoUI
    using StaticArrays
end

# ╔═╡ 0749dd0c-e2ee-4163-b737-2c992966fb0c
setprecision(Arb, 128)

# ╔═╡ 0a71b062-e8d6-4b38-8b15-4b84078ca133
TableOfContents()

# ╔═╡ 031d4b9b-2332-42c1-8adc-a424053bfdac
md"""
## Setup parameters
"""

# ╔═╡ dc4b2b34-366d-4221-a4a2-1c83edf4a823
p = gl_params(1, Arb(1.0), 2.3, 0.0, 0.0)

# ╔═╡ 4d67c703-134d-4d98-81f5-a10df3735342
κ, μ = Arb(0.49323), Arb(0.78308)

# ╔═╡ dd714129-d38e-406a-aa33-d58c3988f66f
ξ₁ = Arb(10)

# ╔═╡ d4a0c14d-d6eb-4396-8580-4447c4416862
md"""
## Constant for $U$
"""

# ╔═╡ d1129f5c-440a-4543-b704-8d889055d611
md"""
The constant $C_U$ (depending on $z_1$) should satisfy

$$|U(a, b, z)| \leq C_U |z^{-a}|$$

for $|z| \geq |z_1|$ and with $|\mathrm{Im}(z)| > |b - a|$.
"""

# ╔═╡ 27db6f6b-68eb-48a0-9f23-f176f0ff798e
CU = let
    d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

    a = Acb(1 / σ, ω / κ) / 2
    b = Acb(d // 2)
    z₁ = Acb(0, -κ) / Acb(1, -ϵ) * ξ₁^2 / 2
    GinzburgLandauSelfSimilarSingular.C_hypgeom_u(a, b, z₁, 5)
end

# ╔═╡ f484d5f3-b2de-4dad-b637-ca98aa81bf1a
let
    pl = plot()

    ξs = range(ξ₁, 100, 100)

    a = Acb(1 / p.σ, p.ω / κ) / 2
    b = Acb(p.d // 2)
    zs = map(ξ -> Acb(0, -κ) / Acb(1, -p.ϵ) * ξ^2 / 2, ξs)

    f = z -> begin
        GinzburgLandauSelfSimilarSingular.hypgeom_u(a, b, z) |> abs
    end

    fs = f.(zs)
    hs = (z -> CU * abs(z^-a)).(zs)

    plot!(
        pl,
        ξs,
        abs.(GinzburgLandauSelfSimilarSingular.hypgeom_u.(a, b, zs) .* zs .^ a),
        label = L"|U(a, b, z)|",
    )

    plot!(pl, ξs, CU * abs.(zs .^ -a .* zs .^ a), label = L"C_U |z^{-a}|")

    pl
end

# ╔═╡ 6fafb070-3b7d-4b61-9d89-e96e92733768
md"""
## Constants for $P$ and $E$
"""

# ╔═╡ 02dd25a6-3ef7-4d96-8451-2fb6743b293d
md"""
The constant $C_P$ should satisfy

$$|P(\xi; p, \kappa)| \leq C_P \xi^{-1 / \sigma}$$

for $\xi \geq 1$.
"""

# ╔═╡ 31cbaa43-e6f8-4b6b-82b0-ea9e46e64cec
CP = GinzburgLandauSelfSimilarSingular.C_P(κ, p, ξ₁)

# ╔═╡ 8668f5ee-7426-40e6-a4a9-dee7df34de44
let
    pl = plot()

    plot!(
        pl,
        range(Arb(ξ₁), 100, 1000),
        ξ -> abs(P(ξ, (p, κ))) * ξ^(1 / p.σ),
        label = L"|P(\xi; p, \kappa)| \xi^{1 / \sigma}",
    )

    hline!(pl, [CP], label = L"C_P")

    pl
end

# ╔═╡ 047d0e9e-9568-4568-a420-e055cb58734d
md"""
The constant $C_E$ should satisfy

$$|E(\xi; p, \kappa)| \leq C_E e^{\mathrm{Re}(z)} \xi^{-d + 1 / \sigma}$$

for $\xi \geq 1$ and with $z = -\frac{i \kappa}{2(1 - i \epsilon)}\xi^2$. Note that $\mathrm{Re}(z) = \frac{\kappa \epsilon}{2(1 + \epsilon^2)}\xi^2$.
"""

# ╔═╡ ee670e0c-69cd-4343-a7f3-8979bd864d6f
CE = GinzburgLandauSelfSimilarSingular.C_E(κ, p, ξ₁)

# ╔═╡ 259e929b-bd8c-46a9-8425-cf67f2d1fd30
let E = GinzburgLandauSelfSimilarSingular.E

    c = κ * p.ϵ / 2(1 + p.ϵ^2)

    pl = plot()

    plot!(
        pl,
        range(ξ₁, 20, 1000),
        ξ -> abs(E(ξ, (p, κ))) * ξ^(p.d - 1 / p.σ) * exp(-c * ξ^2),
        label = L"|E(\xi; p, \kappa)| \xi^{d - 1 / \sigma} e^{-\mathrm{Re}(z)}",
    )

    hline!(pl, [CE], label = L"C_E")
    #ylims!(pl, 0, NaN)

    pl
end

# ╔═╡ 550b34af-8c21-44f0-9ee5-f147a64964f8
md"""
## Constant for $K$
"""

# ╔═╡ 47e7d91e-2c10-44cb-a34e-054686f29b8e
md"""
The constant $C_K$ should satisfy

$$|K(\xi, \eta; p, \kappa)| \leq C_K \xi^{-1 / \sigma} \nu^{1 / \sigma - 1}$$

for $1 \leq \eta \leq \xi$ and

$$|K(\xi, \eta; p, \kappa)| \leq C_K \xi^{-d + 1 / \sigma} \nu^{-1 - 1 / \sigma + d}$$

for $1 \leq \xi \leq \eta$.
"""

# ╔═╡ 65d1f85d-dfdd-40ba-b8c9-b86b92728209
CK = GinzburgLandauSelfSimilarSingular.C_K(κ, p, ξ₁)

# ╔═╡ cc4e4c23-49dd-4980-b65f-4f4fed7f9d5d
let K = GinzburgLandauSelfSimilarSingular.K

    pl = plot()

    factor = (ξ, η) -> if η <= ξ
        ξ^(-1 / p.σ) * η^(1 / p.σ - 1)
    else
        ξ^(-p.d + 1 / p.σ) * η^(-1 - 1 / p.σ + p.d)
    end

    for η in [10, 15, 20, 50]
        plot!(
            pl,
            range(ξ₁, 50, 250),
            ξ -> abs(K(ξ, Arb(η), (p, κ))) / factor(ξ, Arb(η)),
            label = "η = $η",
        )
    end

    hline!(pl, [CK], label = L"C_E")

    pl
end

# ╔═╡ 50cf880b-f997-4942-9c75-ba0a538dadaf
md"""
## Constants for $T$
"""

# ╔═╡ 937364b9-e284-4ab3-96a6-09fdb12172a3
md"""
The constants $C_{T,1}$ and $C_{T,2}$ should satisfy

$$\|T(u; p, \kappa)\|_{v} \leq C_{T,1}|\gamma|\xi_1^{-v} + C_{T,2}\xi_1^{-2 + 2\sigma v}\|u\|_{v}^{2\sigma + 1}$$

where $\xi_1 \geq 1$, $\|u\|_v$ is the norm given by

$$\|u\|_{v} = \sup_{\xi \geq \xi_1} \xi^{1 / \sigma - v}|u(\xi)|$$

and $T$ is the operator given by

$$T(u; p, \kappa) = \gamma P(\xi; p, \kappa) - \int_{\xi_1}^\infty (1 + i\delta)K(\xi, \eta; p, \kappa)|u(\eta)|^{2\sigma}u(\eta)\ d\eta.$$

"""

# ╔═╡ 59cabc8a-5f85-4b39-a35e-ae145bb4a513
let
    pl = plot(xlabel = "v")

    f = v -> try
        GinzburgLandauSelfSimilarSingular.C_T1(v, κ, p, ξ₁)[2]
    catch
        GinzburgLandauSelfSimilarSingular.indeterminate(v)
    end

    plot!(pl, range(Arb(0), 0.5, 100), f, label = L"C_{T,2}")

    pl
end

# ╔═╡ 81f77a06-4f41-4604-bdcf-d4bc963cc03d
md"""
The constant $C_{T,3}$ should satisfy 

$$\|T(u_1; p, \kappa) - T(u_2; p, \kappa)\|_v \leq C_{T,3}\xi_1^{-2 + 2\sigma v}\|u_1 - u_2\|_v(\|u_1\|_v^{2\sigma} + \|u_2\|_v^{2\sigma}).$$
"""

# ╔═╡ 0b80348c-39cb-457f-b234-2063373c150e
let
    pl = plot(xlabel = "v")

    f = v -> try
        GinzburgLandauSelfSimilarSingular.C_T2(v, κ, p, ξ₁)
    catch
        GinzburgLandauSelfSimilarSingular.indeterminate(v)
    end

    plot!(pl, range(Arb(0), 0.5, 100), f, label = L"C_{T,3}")

    pl
end

# ╔═╡ 807674fb-0551-4035-81f2-89f375cf2059
md"""
## Valid values for $r_1$ and $\rho$
If we have that

$$C_{T,1}r_1\xi_1^{-v} + C_{T,2}\xi_1^{-2 + 2\sigma v}\rho^{2\sigma + 1} \leq \rho$$

and

$$2C_{T,3}\rho^{2\sigma}\xi_1^{-2 + 2\sigma v} < 1$$

then the operator $T$ is a contraction of the ball

$$B_\rho = \{u : \|u\|_v \leq \rho\}$$

for $|\gamma| \leq r_1$.
"""

# ╔═╡ 64b553ce-9a79-4e55-bdd8-8b7e55ecd2c1
md"""
For given $v$, $\xi_1$, $\rho$ and $r_1$ we can check if the inequalities are satisfied.
"""

# ╔═╡ e5709fbd-13e2-4bec-8498-e27086b237ce
function is_feasible(v, ξ₁, ρ, r₁; CT₁ = nothing, CT₂ = nothing, CT₃ = nothing)

    if isnothing(CT₁)
        CT₁, CT₂, CT₃ = try
            CT₁, CT₂ = GinzburgLandauSelfSimilarSingular.C_T1(v, κ, p, ξ₁)
            CT₃ = GinzburgLandauSelfSimilarSingular.C_T2(v, κ, p, ξ₁)
            CT₁, CT₂, CT₃
        catch
            return false
        end
    end

    return CT₁ * r₁ * ξ₁^(-v) + CT₂ * ξ₁^(-2 + 2p.σ * v) * ρ^(2p.σ + 1) <= ρ &&
           2CT₃ * ρ^(2p.σ) * ξ₁^(-2 + 2p.σ * v) < 1
end

# ╔═╡ 409d99c7-8374-4a3d-998d-709a7b87e17e
md"""
For fixed $v$ and $\xi_1$ we can then plot the set of feasible $\rho$ and $r_1$.
"""

# ╔═╡ bc70b51e-ddd3-4f47-b0eb-0b4cf66d94f5
@bind v_slider Slider(range(0, 0.4, 100), default = 0.1, show_value = true)

# ╔═╡ 9c3a8515-acdf-42f8-8539-34f7195c3731
@bind ξ₁_slider Slider(range(1, 100, 100), default = 30, show_value = true)

# ╔═╡ eb1f56dc-c92e-4a6e-9c42-0466a8f0019b
let v = Arb(v_slider), ξ₁ = Arb(ξ₁_slider)
    pl = plot(xlabel = L"\rho", ylabel = L"r_1", colorbar = :none)

    CT₁, CT₂ = GinzburgLandauSelfSimilarSingular.C_T1(v, κ, p, ξ₁)
    CT₃ = GinzburgLandauSelfSimilarSingular.C_T2(v, κ, p, ξ₁)

    ρ_max = (ξ₁^(2 - 2p.σ * v) / 2CT₃)^(1 / 2p.σ)

    ρs = range(Arb(0), 1.05ρ_max, 100)
    r₁s = range(Arb(0), 50, 100)

    res = map((ρ, r₁) for r₁ in r₁s, ρ in ρs) do (ρ, r₁)
        is_feasible(v, ξ₁, ρ, r₁; CT₁, CT₂, CT₃)
    end

    heatmap!(pl, ρs, r₁s, res)

    pl
end

# ╔═╡ fecf0b3a-bbc0-40e7-a17a-450b43643494
let v = Arb(v_slider), ξ₁ = Arb(ξ₁_slider)
    pl = plot(xlabel = L"\rho", ylabel = L"r_1", colorbar = :none)

    CT₁, CT₂ = GinzburgLandauSelfSimilarSingular.C_T1(v, κ, p, ξ₁)
    CT₃ = GinzburgLandauSelfSimilarSingular.C_T2(v, κ, p, ξ₁)

    ρ_max = (ξ₁^(2 - 2p.σ * v) / 2CT₃)^(1 / 2p.σ)

    ρs = range(Arb(0), ρ_max, 100)
    r₁s = map(ρs) do ρ
        (ρ - CT₂ * ξ₁^(-2 + 2p.σ * v) * ρ^(2p.σ + 1)) * ξ₁^v / CT₁
    end

    plot!(ρs, r₁s, fillrange = [0])

    pl
end

# ╔═╡ 808dbac2-fbd7-41c0-8e6c-22171ee81ae2
md"""
## Bounds for fix point
The constants $C_{\infty,1}$ and $C_{\infty,2}$ should satisfy that for the fix point $u_\infty$ of the operator $T$ we have

$$u_\infty(\xi) = \gamma P(\xi) + P(\xi)f(\xi) + E(\xi)g(\xi)$$

where 

$$|f(\xi)| \leq C_{\infty,1}(\xi_1^{2\sigma v + v - 2} - \xi^{2\sigma v + v - 2})$$ 

and 

$$|g(\xi)| \leq C_{\infty,2}e^{-\mathrm{Re}(c)\xi^2}\xi^{-2 / \sigma + 2\sigma v + v + d - 2}$$

for $\xi \geq \xi_1$ and $|\gamma| \leq r_1$. In particular for $\xi = \xi_1$ we then get

$$u_\infty(\xi) = \gamma P(\xi) + E(\xi)g(\xi).$$
"""

# ╔═╡ cd9dcd58-18cd-41fe-ad12-21cd6a630d05
md"""
The constants $C_{\infty,3}$ and $C_{\infty,4}$ together with the two above constants should satisfy

$$u_\infty'(\xi) = \gamma P'(\xi) + P'(\xi)f(\xi) + P(\xi)f'(\xi) + E'(\xi)g(\xi) + E(\xi)g'(\xi)$$

where

$$|f'(\xi)| \leq C_{\infty,3}\xi^{-2\sigma v + v - 3}$$

and

$$|g'(\xi)| \leq C_{\infty,4}e^{-\mathrm{Re}(c)\xi^2}\xi^{-2 / \sigma + 2\sigma v + v + d - 3}$$
"""

# ╔═╡ d39e45eb-cb96-478c-9912-26a2fdfc241e
@bind v_slider_2 Slider(range(0, 0.4, 100), default = 0.1, show_value = true)

# ╔═╡ 2263b014-818d-46fa-a37f-886a036f1a90
@bind ξ₁_slider_2 Slider(range(1, 100, 100), default = 30, show_value = true)

# ╔═╡ 5f3a4c77-5caa-4016-a6f4-cad67139978f
@bind r₁_slider_2 Slider(range(0, 10, 100), default = 1, show_value = true)

# ╔═╡ 51a52721-0e73-4336-b730-83eb05984590
C∞₁, C∞₂ = let v = Arb(v_slider_2), ξ₁ = Arb(ξ₁_slider_2), r₁ = Arb(r₁_slider_2)
    GinzburgLandauSelfSimilarSingular.C_fix_point(r₁, v, κ, p, ξ₁)
end

# ╔═╡ cbb68274-7a11-4555-b19e-ce9062610fd7
C∞₃, C∞₄ = let v = Arb(v_slider_2), ξ₁ = Arb(ξ₁_slider_2), r₁ = Arb(r₁_slider_2)
    GinzburgLandauSelfSimilarSingular.C_fix_point_derivative(r₁, v, κ, p, ξ₁)
end

# ╔═╡ b110a15b-7d8c-4c15-b9a6-692e20f286ae
md"""
### Bounds for explicit $\gamma$
"""

# ╔═╡ bb9d51e5-c460-4364-92f0-9c7205e28933
let v = Arb(v_slider_2), ξ₁ = Arb(ξ₁_slider_2)#, r₁ = Arb(r₁_slider_2)
    γ = Acb(1.7920178776656215, 1.1043012090074602)
    r₁ = 1.01abs(γ)

    C∞₁, C∞₂ = GinzburgLandauSelfSimilarSingular.C_fix_point(r₁, v, κ, p, ξ₁)
    C∞₃, C∞₄ = GinzburgLandauSelfSimilarSingular.C_fix_point_derivative(r₁, v, κ, p, ξ₁)
    rc = κ * p.ϵ / (1 + p.ϵ)^2

    E = GinzburgLandauSelfSimilarSingular.E
    E_dξ = GinzburgLandauSelfSimilarSingular.E_dξ

    ξs = range(ξ₁, ξ₁ + 10, 100)
    us = map(ξs) do ξ
        f_abs_bound = C∞₁ * (ξ₁^(2p.σ * v + v - 2) - ξ^(2p.σ * v + v - 2))
        g_abs_bound = C∞₂ * exp(-rc * ξ^2) * ξ^(-2 / p.σ + 2p.σ * v + v + p.d - 2)

        f_ball = add_error(Acb(), f_abs_bound)
        g_ball = add_error(Acb(), g_abs_bound)

        γ * P(ξ, (p, κ)) + P(ξ, (p, κ)) * f_ball + E(ξ, (p, κ)) * g_ball
    end
    us = abs.(us)

    dus = map(ξs) do ξ
        f_abs_bound = C∞₁ * (ξ₁^(2p.σ * v + v - 2) - ξ^(2p.σ * v + v - 2))
        g_abs_bound = C∞₂ * exp(-rc * ξ^2) * ξ^(-2 / p.σ + 2p.σ * v + v + p.d - 2)
        df_abs_bound = C∞₃ * ξ^(2p.σ * v + v - 3)
        dg_abs_bound = C∞₄ * exp(-rc * ξ^2) * ξ^(-2 / p.σ + 2p.σ * v + v + p.d - 3)

        f_ball = add_error(Acb(), f_abs_bound)
        g_ball = add_error(Acb(), g_abs_bound)
        df_ball = add_error(Acb(), df_abs_bound)
        dg_ball = add_error(Acb(), dg_abs_bound)

        γ * P(ξ, (p, κ)) +
        P_dξ(ξ, (p, κ)) * f_ball +
        P(ξ, (p, κ)) * df_ball +
        E_dξ(ξ, (p, κ)) * g_ball +
        E(ξ, (p, κ)) * dg_ball
    end
    dus = abs.(dus)

    pl = plot(xlabel = L"\xi", ylabel = "Error")

    plot!(
        pl,
        ξs,
        [radius.(Arb, us) radius.(Arb, dus)],
        label = [L"u_\infty(\xi)" L"u_\infty'(\xi)"],
    )
end

# ╔═╡ Cell order:
# ╠═8d7322fe-e018-11ed-1579-7daddd6ce28f
# ╠═0749dd0c-e2ee-4163-b737-2c992966fb0c
# ╠═0a71b062-e8d6-4b38-8b15-4b84078ca133
# ╟─031d4b9b-2332-42c1-8adc-a424053bfdac
# ╠═dc4b2b34-366d-4221-a4a2-1c83edf4a823
# ╠═4d67c703-134d-4d98-81f5-a10df3735342
# ╠═dd714129-d38e-406a-aa33-d58c3988f66f
# ╟─d4a0c14d-d6eb-4396-8580-4447c4416862
# ╟─d1129f5c-440a-4543-b704-8d889055d611
# ╠═27db6f6b-68eb-48a0-9f23-f176f0ff798e
# ╟─f484d5f3-b2de-4dad-b637-ca98aa81bf1a
# ╟─6fafb070-3b7d-4b61-9d89-e96e92733768
# ╟─02dd25a6-3ef7-4d96-8451-2fb6743b293d
# ╠═31cbaa43-e6f8-4b6b-82b0-ea9e46e64cec
# ╟─8668f5ee-7426-40e6-a4a9-dee7df34de44
# ╟─047d0e9e-9568-4568-a420-e055cb58734d
# ╠═ee670e0c-69cd-4343-a7f3-8979bd864d6f
# ╟─259e929b-bd8c-46a9-8425-cf67f2d1fd30
# ╟─550b34af-8c21-44f0-9ee5-f147a64964f8
# ╟─47e7d91e-2c10-44cb-a34e-054686f29b8e
# ╠═65d1f85d-dfdd-40ba-b8c9-b86b92728209
# ╟─cc4e4c23-49dd-4980-b65f-4f4fed7f9d5d
# ╟─50cf880b-f997-4942-9c75-ba0a538dadaf
# ╟─937364b9-e284-4ab3-96a6-09fdb12172a3
# ╟─59cabc8a-5f85-4b39-a35e-ae145bb4a513
# ╟─81f77a06-4f41-4604-bdcf-d4bc963cc03d
# ╟─0b80348c-39cb-457f-b234-2063373c150e
# ╟─807674fb-0551-4035-81f2-89f375cf2059
# ╟─64b553ce-9a79-4e55-bdd8-8b7e55ecd2c1
# ╠═e5709fbd-13e2-4bec-8498-e27086b237ce
# ╟─409d99c7-8374-4a3d-998d-709a7b87e17e
# ╠═bc70b51e-ddd3-4f47-b0eb-0b4cf66d94f5
# ╠═9c3a8515-acdf-42f8-8539-34f7195c3731
# ╟─eb1f56dc-c92e-4a6e-9c42-0466a8f0019b
# ╟─fecf0b3a-bbc0-40e7-a17a-450b43643494
# ╟─808dbac2-fbd7-41c0-8e6c-22171ee81ae2
# ╟─cd9dcd58-18cd-41fe-ad12-21cd6a630d05
# ╠═d39e45eb-cb96-478c-9912-26a2fdfc241e
# ╠═2263b014-818d-46fa-a37f-886a036f1a90
# ╠═5f3a4c77-5caa-4016-a6f4-cad67139978f
# ╠═51a52721-0e73-4336-b730-83eb05984590
# ╠═cbb68274-7a11-4555-b19e-ce9062610fd7
# ╟─b110a15b-7d8c-4c15-b9a6-692e20f286ae
# ╠═bb9d51e5-c460-4364-92f0-9c7205e28933
