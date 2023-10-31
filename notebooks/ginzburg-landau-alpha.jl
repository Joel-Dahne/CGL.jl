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

# ╔═╡ 60dfdc80-ba5e-11ed-16f4-ade6c584dc7b
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using Arblib
    using DifferentialEquations
    using CGL
    using NLsolve
    using Plots
    using PlutoUI
end

# ╔═╡ 374e1504-4ea4-4c24-be24-ea6a85b0a11e
setprecision(Arb, 512)

# ╔═╡ d3d1fa98-770d-4f14-bfa4-ca0a28a35563
md"""
# The interval $(0, ξ₀)$
"""

# ╔═╡ f57069eb-ece9-4744-aa10-31edfc25e63b
md"""
We want to solve the ODE

$$a'' + ϵb'' + \frac{d - 1}{ξ}(a' + ϵb') - κξb' - \frac{κ}{σ}b - ωa + (a^2 + b^2)^σ a - δ(a^2 + b^2)^σ b = 0$$
$$b'' - ϵa'' + \frac{d - 1}{ξ}(b' - ϵa') + κξa' + \frac{κ}{σ}a - ωb + (a^2 + b^2)^σ b + δ(a^2 + b^2)^σ a = 0$$

The first step is to rewrite it as a system of first order ODE:s. For that we introduce the variables $x = a'$ and $y = b'$, we then have

$$a' = x$$

$$b' = y$$

$$x' + ϵy'= -\frac{d - 1}{ξ}(x + ϵy) + κξy + \frac{κ}{σ}b + ωa - (a^2 + b^2)^σ a + δ(a^2 + b^2)^σ b$$

$$-ϵx' + y' = -\frac{d - 1}{ξ}(y - ϵx) - κξx - \frac{κ}{σ}a + ωb - (a^2 + b^2)^σ b - δ(a^2 + b^2)^σ a$$

If we let

$$F₁(a, b, x, y, ξ) = -\frac{d - 1}{ξ}(x + ϵy) + κξy + \frac{κ}{σ}b + ωa - (a^2 + b^2)^σ a + δ(a^2 + b^2)^σ b$$

$$F₂(a, b, x, y, ξ) = -\frac{d - 1}{ξ}(y - ϵx) - κξx - \frac{κ}{σ}a + ωb - (a^2 + b^2)^σ b - δ(a^2 + b^2)^σ a$$

When we can solve for $x'$ and $y'$ to get

$$x' = F₁(a, b, x, y, ξ) - ϵF₂(a, b, x, y, ξ)$$

$$y' = ϵF₁(a, b, x, y, ξ) + F₂(a, b, x, y, ξ)$$
"""

# ╔═╡ 15e38430-ca56-4b5b-90c0-0eddc6725329
function f!(du, u, p, ξ)
    a, b, x, y = u
    if p isa NamedTuple
        (; d, ω, κ, σ, ϵ, δ) = p
    else
        d, ω, κ, σ, ϵ, δ = p
    end

    F1 = κ * ξ * y + κ / σ * b + ω * a - (a^2 + b^2)^σ * a + δ * (a^2 + b^2)^σ * b
    F2 = -κ * ξ * x - κ / σ * a + ω * b - (a^2 + b^2)^σ * b - δ * (a^2 + b^2)^σ * a

    if !(x == y == ξ == 0)
        F1 += -(d - 1) / ξ * (x + ϵ * y)
        F2 += -(d - 1) / ξ * (y - ϵ * x)
    end

    du[1] = x
    du[2] = y
    du[3] = F1 - ϵ * F2
    du[4] = ϵ * F1 + F2

    du
end

# ╔═╡ 87cce4db-37c5-43f3-99b6-70467486633f
md"""
## Parameters
"""

# ╔═╡ a79b3baa-ffa8-440a-ab2d-1d7c2d068b71
tspan = (0.0, 30.0)

# ╔═╡ 8042cafa-5e7c-49bb-b646-d0f2a92bf710
md"""
These are the parameters used in Figures 2 to 7
"""

# ╔═╡ dc647047-94e6-4670-a98e-fc07b141f0ea
parameters1 = (d = 1.0, ω = 1.0, κ = 0.49323, σ = 2.3, μ = 0.78308, ϵ = 0.0, δ = 0.0)

# ╔═╡ d0b548ed-ca9f-4d0a-9eb0-e3a6fdfcc689
parameters2 = (d = 1.0, ω = 3.07959, κ = 1.51894, σ = 2.3, μ = 1.0, ϵ = 0.0, δ = 0.0)

# ╔═╡ 268e303d-b41f-4836-9a83-ae1cb0c94eec
parameters3 = (d = 1.0, ω = 1.0, κ = 0.26678, σ = 1.0, μ = 0.78308, ϵ = 0.0, δ = 0.0)

# ╔═╡ 155768b4-5daf-457d-b0d7-c90e37dd52cd
parameters4 = (d = 1.0, ω = 1.76651, κ = 0.47127, σ = 2.3, μ = 1.0, ϵ = 0.0, δ = 0.0)

# ╔═╡ 23dcf903-caf8-47f5-a7ff-38b5df47c6bc
parameters5 = (d = 3.0, ω = 1.0, κ = 0.32091, σ = 1.0, μ = 0.83559, ϵ = 0.0, δ = 0.0)

# ╔═╡ fd138fb8-b459-4888-b12d-89a838c710ff
parameters6 = (d = 3.0, ω = 1.41727, κ = 0.45535, σ = 2.3, μ = 1.0, ϵ = 0.0, δ = 0.0)

# ╔═╡ 9c8042a5-0a86-457d-9de7-f2e059c5764c
md"""
## Reproduce figures
We can now reproduce Figure 2 to 7 from the thesis. 
"""

# ╔═╡ 47b35421-7db1-4738-a20d-d79b4801cd36
sol1 = let p = parameters1
    prob = ODEProblem(f!, [p.μ, 0.0, 0.0, 0.0], tspan, p)
    solve(prob, abstol = 1e-9, reltol = 1e-9)
end;

# ╔═╡ 7394e042-228a-4915-abd4-dc00135b2f8a
sol2 = let p = parameters2
    prob = ODEProblem(f!, [p.μ, 0.0, 0.0, 0.0], tspan, p)
    solve(prob, abstol = 1e-9, reltol = 1e-9)
end;

# ╔═╡ 07763f0c-23eb-476c-a818-3767d72df5fc
sol3 = let p = parameters3
    prob = ODEProblem(f!, [p.μ, 0.0, 0.0, 0.0], tspan, p)
    solve(prob, abstol = 1e-9, reltol = 1e-9)
end;

# ╔═╡ e38e45bd-cb13-4fbb-a4e0-e80aff2c48c2
sol4 = let p = parameters4
    prob = ODEProblem(f!, [p.μ, 0.0, 0.0, 0.0], tspan, p)
    solve(prob, abstol = 1e-9, reltol = 1e-9)
end;

# ╔═╡ b2a0efc4-67ad-4fe7-8813-93d83bb79805
sol5 = let p = parameters5
    prob = ODEProblem(f!, [p.μ, 0.0, 0.0, 0.0], tspan, p)
    solve(prob, abstol = 1e-9, reltol = 1e-9)
end;

# ╔═╡ 75640143-513d-42bc-af58-55858cc4086a
sol6 = let p = parameters6
    prob = ODEProblem(f!, [p.μ, 0.0, 0.0, 0.0], tspan, p)
    solve(prob, abstol = 1e-9, reltol = 1e-9)
end;

# ╔═╡ b4299bd3-6a9c-4216-83c3-34495feed058
let sol = sol1
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])
    plot(pl1, pl2)
end

# ╔═╡ 03f2f89c-8a8e-45d1-9dd5-05d12ad6e3a6
let sol = sol2
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])
    plot(pl1, pl2)
end

# ╔═╡ 27685362-25d8-4f4a-a5bf-7787d9760aaf
let sol = sol3
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])
    plot(pl1, pl2, title = "Not like in thesis!")
end

# ╔═╡ 85fc701b-c12e-4d0a-a404-2c214e23d979
let sol = sol4
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])
    plot(pl1, pl2)
end

# ╔═╡ 150d8a06-7c05-40d5-a03d-5ebff8b04558
let sol = sol5
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])
    plot(pl1, pl2)
end

# ╔═╡ 89c7946b-3d56-41d0-a9dc-170061f246be
let sol = sol6
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])
    plot(pl1, pl2, title = "This one seems wrong!")
end

# ╔═╡ 6b764bdc-cc98-4411-a021-0fea807b7a11
md"""
## Poincaré-Miranda Theorem
"""

# ╔═╡ ff26580d-c378-46b9-9058-7d03ae67b273
md"""
Define the box we want to apply the Poincaré-Miranda Theorem on.
"""

# ╔═╡ b04575eb-25a5-4839-9ef6-213429e7485a
begin
    n = 10
    ω_endpoints = (1.0, 1.5)
    κ_endpoints = (0.5, 1.0)
    σ_endpoints = (2.0, 2.3)
    μ_endpoints = (0.5, 1.0)
    ωs = range(ω_endpoints..., n)
    κs = range(κ_endpoints..., n)
    σs = range(σ_endpoints..., n)
    μs = range(μ_endpoints..., n)
end

# ╔═╡ e8f0a6e2-81c9-44c9-8b88-6bb44c42d1d6
md"""
Compute values for the different faces of the box.
"""

# ╔═╡ 4780850f-e1d0-434f-a25a-e4501a5fef84
ω_faces = let
    faces = (Array{Vector{Float64}}(undef, n, n, n), Array{Vector{Float64}}(undef, n, n, n))
    for (ω, face) in zip(ω_endpoints, faces)
        for (i, κ) in enumerate(κs)
            for (j, σ) in enumerate(σs)
                for (k, μ) in enumerate(μs)
                    p = [1.0, ω, κ, σ, 0.0, 0.0]
                    u0 = [μ, 0.0, 0.0, 0.0]
                    tspan = (0.0, 30.0)
                    prob = ODEProblem(f!, u0, tspan, p)
                    sol = solve(prob)
                    face[i, j, k] = sol[end]
                end
            end
        end
    end
    faces
end

# ╔═╡ 1828428f-c5e8-4af5-84e2-42a504eaa386
κ_faces = let
    faces = (Array{Vector{Float64}}(undef, n, n, n), Array{Vector{Float64}}(undef, n, n, n))
    for (κ, face) in zip(κ_endpoints, faces)
        for (i, ω) in enumerate(ωs)
            for (j, σ) in enumerate(σs)
                for (k, μ) in enumerate(κs)
                    p = [1.0, ω, κ, σ, 0.0, 0.0]
                    u0 = [μ, 0.0, 0.0, 0.0]
                    tspan = (0.0, 30.0)
                    prob = ODEProblem(f!, u0, tspan, p)
                    sol = solve(prob)
                    face[i, j, k] = sol[end]
                end
            end
        end
    end
    faces
end

# ╔═╡ eebe6cd5-2a1b-484b-acf8-ea117ee6ef3c
σ_faces = let
    faces = (Array{Vector{Float64}}(undef, n, n, n), Array{Vector{Float64}}(undef, n, n, n))
    for (σ, face) in zip(σ_endpoints, faces)
        for (i, ω) in enumerate(ωs)
            for (j, κ) in enumerate(κs)
                for (k, μ) in enumerate(μs)
                    p = [1.0, ω, κ, σ, 0.0, 0.0]
                    u0 = [μ, 0.0, 0.0, 0.0]
                    tspan = (0.0, 30.0)
                    prob = ODEProblem(f!, u0, tspan, p)
                    sol = solve(prob)
                    face[i, j, k] = sol[end]
                end
            end
        end
    end
    faces
end

# ╔═╡ e78dc117-66a7-4ea0-a487-0d8591d60107
μ_faces = let
    faces = (Array{Vector{Float64}}(undef, n, n, n), Array{Vector{Float64}}(undef, n, n, n))
    for (μ, face) in zip(μ_endpoints, faces)
        for (i, ω) in enumerate(ωs)
            for (j, κ) in enumerate(κs)
                for (k, σ) in enumerate(σs)
                    p = [1.0, ω, κ, σ, 0.0, 0.0]
                    u0 = [μ, 0.0, 0.0, 0.0]
                    tspan = (0.0, 30.0)
                    prob = ODEProblem(f!, u0, tspan, p)
                    sol = solve(prob)
                    face[i, j, k] = sol[end]
                end
            end
        end
    end
    faces
end

# ╔═╡ e34216ff-716c-42f1-9bb2-f66117da6446
md"""
In the end we will be interested in the sign of the k-th component at the k-th face.
"""

# ╔═╡ 15315635-2762-475f-8e5f-05361f51a483
ω_faces_sign_a = Tuple(sign.(getindex.(face, 1)) for face in ω_faces)

# ╔═╡ 9e936fa9-8953-4104-becf-aa4ffc14571d
κ_faces_sign_b = Tuple(sign.(getindex.(face, 2)) for face in κ_faces)

# ╔═╡ 1eb3e568-37f6-43a4-86ae-686d0ce84bc4
σ_faces_sign_x = Tuple(sign.(getindex.(face, 3)) for face in σ_faces)

# ╔═╡ 13dc0e82-1a67-44a8-9060-37dc986ba518
μ_faces_sign_y = Tuple(sign.(getindex.(face, 4)) for face in μ_faces)

# ╔═╡ 6a0dce19-a6ff-4e1f-a332-bf558d0e37be
md"""
We can now check the requirement of the theorem. Note that for the values we currently compute it will not be satisfied.
"""

# ╔═╡ a9e9a92f-dc53-492d-bc86-9753356c84ec
theorem_ω_requirement =
    all(<=(0), ω_faces_sign_a[1]) && all(>=(0), ω_faces_sign_a[2]) ||
    all(>=(0), ω_faces_sign_a[1]) && all(<=(0), ω_faces_sign_a[2])

# ╔═╡ 14111c0e-89c0-41cd-92bd-cf9c080a2913
theorem_κ_requirement =
    all(<=(0), κ_faces_sign_b[1]) && all(>=(0), κ_faces_sign_b[2]) ||
    all(>=(0), κ_faces_sign_b[1]) && all(<=(0), κ_faces_sign_b[2])

# ╔═╡ cc684d4f-0feb-4150-a446-db664c56f678
theorem_σ_requirement =
    all(<=(0), σ_faces_sign_x[1]) && all(>=(0), σ_faces_sign_x[2]) ||
    all(>=(0), σ_faces_sign_x[1]) && all(<=(0), σ_faces_sign_x[2])

# ╔═╡ e0a50a81-6312-4607-89f3-5f44b8b3d568
theorem_μ_requirement =
    all(<=(0), μ_faces_sign_y[1]) && all(>=(0), μ_faces_sign_y[2]) ||
    all(>=(0), μ_faces_sign_y[1]) && all(<=(0), μ_faces_sign_y[2])

# ╔═╡ 4ed0cb0c-485b-4d2f-b028-3aff709f37d9
md"""
Try making some plots for visualization, this is difficualy due to the high dimensionionality.
"""

# ╔═╡ 8125af3e-ee25-4b84-9279-5d67f95d1772
scatter3d(κs, ωs, getindex.(ω_faces[1][:, :, 1], 1))

# ╔═╡ c09d1ce4-9012-4cea-94a4-06a255b93033
getindex.(ω_faces[1][:, :, 1], 1)

# ╔═╡ 712552b8-73e0-4245-b12c-3ac63ed13e14
md"""
## Taylor seriers
"""

# ╔═╡ ce107d7b-2c65-4b5c-9113-d1e293c2a846
md"""
## Taylor series at zero
First step is to be able to compute Taylor series at $ξ = 0$. This uses Equation 3.4 in the thesis.
"""

# ╔═╡ 322e24da-bdfd-4298-92a6-651043a3c03f
function Q_expansion_zero((; d, ω, κ, σ, μ, ϵ, δ); degree = 5)
    a = ArbSeries([μ, 0]; degree)
    b = ArbSeries([0, 0]; degree)

    for n = 0:degree-2
        u1 = (a^2 + b^2)^σ * a
        u2 = (a^2 + b^2)^σ * b
        F = κ * n * b[n] + κ / σ * b[n] + ω * a[n] - u1[n] + δ * u2[n]
        G = -κ * n * a[n] - κ / σ * a[n] + ω * b[n] - u2[n] - δ * u1[n]
        a[n+2], b[n+2] = [1 ϵ; -ϵ 1] \ [F; G] / ((n + 2) * (n + d))
    end

    return a, b
end

# ╔═╡ 3035f690-7aff-4192-b8f2-a1e78094047c
a1, b1 = Q_expansion_zero(parameters1, degree = 10)

# ╔═╡ b39e0201-6458-4f87-bbe4-f41a5e40adb0
a2, b2 = Q_expansion_zero(parameters2, degree = 50)

# ╔═╡ 568061ef-fa5d-45d7-9c59-1a0592099dbd
a3, b3 = Q_expansion_zero(parameters3, degree = 20)

# ╔═╡ 3da81649-a810-45ac-b520-aa7ada1fb5db
a4, b4 = Q_expansion_zero(parameters4, degree = 10)

# ╔═╡ 9b18d3eb-5679-428c-98bc-41d12936541f
a5, b5 = Q_expansion_zero(parameters5, degree = 10)

# ╔═╡ 453ee689-ba8c-475c-814d-2bc120e35cc7
a6, b6 = Q_expansion_zero(parameters6, degree = 10)

# ╔═╡ 2ca1a635-5cf1-4aa7-a0ed-692cd87369f1
let a = a1, b = b1, sol = sol1
    ts = range(0, 1, 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    plot!(pl1, ts, abs.(a.(ts) + im * b.(ts)), label = "", linestyle = :dash)
    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    plot!(pl2, ts, [a.(ts) b.(ts)], label = "", color = [:red :blue], linestyle = :dash)

    plot(pl1, pl2)
end

# ╔═╡ 06470b3c-df86-48e3-a69b-9a2c27921c57
let a = a2, b = b2, sol = sol2
    ts = range(0, 1, 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    plot!(pl1, ts, abs.(a.(ts) + im * b.(ts)), label = "", linestyle = :dash)
    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    plot!(pl2, ts, [a.(ts) b.(ts)], label = "", color = [:red :blue], linestyle = :dash)

    plot(pl1, pl2)
end

# ╔═╡ 05a54cee-5829-4cca-94a8-0a408a4c9716
let a = a3, b = b3, sol = sol3
    ts = range(0, 1.0, 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    plot!(pl1, ts, abs.(a.(ts) + im * b.(ts)), label = "", linestyle = :dash)
    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    plot!(pl2, ts, [a.(ts) b.(ts)], label = "", color = [:red :blue], linestyle = :dash)

    plot(pl1, pl2)
end

# ╔═╡ 5bf7172d-84ff-4da1-9e7a-e44d7fed9368
let a = a4, b = b4, sol = sol4
    ts = range(0, 1, 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    plot!(pl1, ts, abs.(a.(ts) + im * b.(ts)), label = "", linestyle = :dash)
    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    plot!(pl2, ts, [a.(ts) b.(ts)], label = "", color = [:red :blue], linestyle = :dash)

    plot(pl1, pl2)
end

# ╔═╡ c9a6df8d-31c4-4e86-88e9-c4c6ba426e3d
let a = a5, b = b5, sol = sol5
    ts = range(0, 1, 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    plot!(pl1, ts, abs.(a.(ts) + im * b.(ts)), label = "", linestyle = :dash)
    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    plot!(pl2, ts, [a.(ts) b.(ts)], label = "", color = [:red :blue], linestyle = :dash)

    plot(pl1, pl2)
end

# ╔═╡ b62d1e45-c91f-4d5e-8e0c-9001fe0481f9
let a = a6, b = b6, sol = sol6
    ts = range(0, 1, 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    plot!(pl1, ts, abs.(a.(ts) + im * b.(ts)), label = "", linestyle = :dash)
    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    plot!(pl2, ts, [a.(ts) b.(ts)], label = "", color = [:red :blue], linestyle = :dash)

    plot(pl1, pl2)
end

# ╔═╡ e4a7ee1f-266b-4553-92b1-6f168f63c9c2
md"""
## Taylor series at a general point
Next step is to compute Taylor series for $ξ > 0$. This uses the equation below Equation 3.9. I believe it should be $n + d$ and not $n + 1$ in the equation though.
"""

# ╔═╡ 3bd6308d-6a33-4a5d-9f03-d50f731f1d3e
function Q_expansion(ξ0, a0, b0, (; d, ω, κ, σ, μ, ϵ, δ); degree = 5)
    a = ArbSeries(a0; degree)
    b = ArbSeries(b0; degree)

    for n = 0:degree-2
        u1 = (a^2 + b^2)^σ * a
        u2 = (a^2 + b^2)^σ * b
        v1 = Arblib.derivative(a) / ArbSeries((ξ0, 1), degree = n)
        v2 = Arblib.derivative(b) / ArbSeries((ξ0, 1), degree = n)
        F =
            (1 - d) * (v1[n] + ϵ * v2[n]) +
            κ * (ξ0 * (n + 1) * b[n+1] + n * b[n]) +
            κ / σ * b[n] +
            ω * a[n] - u1[n] + δ * u2[n]
        G =
            (1 - d) * (v2[n] - ϵ * v1[n]) - κ * (ξ0 * (n + 1) * a[n+1] + n * a[n]) -
            κ / σ * a[n] + ω * b[n] - u2[n] - δ * u1[n]
        a[n+2], b[n+2] = [1 ϵ; -ϵ 1] \ [F; G] / ((n + 2) * (n + 1))
    end

    return a, b
end

# ╔═╡ 0eb0b383-2627-42f2-8038-0ca1241a7735
a11, b11, ξ01 = let ξ0 = 0.5, a = a1, b = b1
    Q_expansion(
        ξ0,
        Arblib.evaluate2(a, ξ0),
        Arblib.evaluate2(b, ξ0),
        parameters1,
        degree = 20,
    )...,
    ξ0
end

# ╔═╡ e8f7d17b-8f6f-43ba-aec7-9fb27a5af8db
a21, b21, ξ02 = let ξ0 = 0.8, a = a2, b = b2
    Q_expansion(
        ξ0,
        Arblib.evaluate2(a, ξ0),
        Arblib.evaluate2(b, ξ0),
        parameters2,
        degree = 10,
    )...,
    ξ0
end

# ╔═╡ 7820cb27-0fb8-44be-a4aa-dd9dbf72b999
let a = a1, b = b1, c = a11, d = b11, ξ0 = ξ01, sol = sol1
    ts = range(0, 2, 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    plot!(pl1, ts, abs.(a.(ts) + im * b.(ts)), label = "", linestyle = :dash)
    plot!(pl1, ts, abs.(c.(ts .- ξ0) + im * d.(ts .- ξ0)), label = "", linestyle = :dot)
    #vline!(pl1, [ξ0])

    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    plot!(pl2, ts, [a.(ts) b.(ts)], label = "", color = [:red :blue], linestyle = :dash)
    plot!(
        pl2,
        ts,
        [c.(ts .- ξ0) d.(ts .- ξ0)],
        label = "",
        color = [:red :blue],
        linestyle = :dot,
    )
    #vline!(pl2, [ξ0])

    plot(pl1, pl2)
end

# ╔═╡ ea5937b1-4ce5-4152-839e-53e8588008dc
let a = a2, b = b2, c = a21, d = b21, ξ0 = ξ02, sol = sol2
    ts = range(0, 1.5, 1000)
    ts1 = range(0, 1.0, 1000)
    ts2 = range(0.5, 1.3, 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    plot!(pl1, ts1, abs.(a.(ts1) + im * b.(ts1)), label = "", linestyle = :dash)
    plot!(pl1, ts2, abs.(c.(ts2 .- ξ0) + im * d.(ts2 .- ξ0)), label = "", linestyle = :dot)
    #vline!(pl1, [ξ0])

    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    plot!(pl2, ts1, [a.(ts1) b.(ts1)], label = "", color = [:red :blue], linestyle = :dash)
    plot!(
        pl2,
        ts2,
        [c.(ts2 .- ξ0) d.(ts2 .- ξ0)],
        label = "",
        color = [:red :blue],
        linestyle = :dot,
    )
    #vline!(pl2, [ξ0])

    plot(pl1, pl2)
end

# ╔═╡ 27e53267-e7aa-4315-b73b-b0a2ad4d0ff1
md"""
## Taylor series at succesive points
Finally we want to iteratitevly use the above methods to go further and further.
"""

# ╔═╡ 4c79b2cc-67ca-412c-acdc-76e9c659b029
function Q_expansions(tspan, Δt, parameters; degree = 10)
    ts = range(tspan..., step = Δt)
    as, bs = similar(ts, ArbSeries), similar(ts, ArbSeries)
    as[1], bs[1] = Q_expansion_zero(parameters; degree)
    for n = 2:length(ts)
        as[n], bs[n] = Q_expansion(
            ts[n],
            Arblib.evaluate2(as[n-1], ts[n] - ts[n-1]),
            Arblib.evaluate2(bs[n-1], ts[n] - ts[n-1]),
            parameters;
            degree,
        )
    end
    return ts, as, bs
end

# ╔═╡ 032f778d-560e-4ea9-91fa-509a5797927b
t0s1, as1, bs1 = Q_expansions((0.0, 30.0), 0.5, parameters1, degree = 20)

# ╔═╡ d267c8ae-dceb-448e-b8b6-33866280cfa5
t0s2, as2, bs2 = Q_expansions((0.0, 15.0), 0.25, parameters2, degree = 20)

# ╔═╡ d7a56d9a-1349-4a0b-bebf-72c9d0afe086
t0s3, as3, bs3 = Q_expansions((0.0, 30.0), 0.25, parameters3, degree = 20)

# ╔═╡ 077462c9-4806-4b22-8cfc-4008d72ff819
t0s4, as4, bs4 = Q_expansions((0.0, 30.0), 0.25, parameters4, degree = 20)

# ╔═╡ dde47573-c285-4e10-9b17-2aae6bd8f125
t0s5, as5, bs5 = Q_expansions((0.0, 30.0), 0.5, parameters5, degree = 20)

# ╔═╡ 381b22aa-87dd-4805-b5b6-d16e2e0e6d8c
t0s6, as6, bs6 = Q_expansions((0.0, 30.0), 0.25, parameters6, degree = 20)

# ╔═╡ 4b3c9483-1483-4335-b101-5a8115157718
let t0s = t0s1, as = as1, bs = bs1, sol = sol1
    Float64[as[end][0], bs[end][0], as[end][1], bs[end][1]] - sol(t0s[end])
end

# ╔═╡ ed7ad5d0-85e3-4994-a073-0c3c78c6cd7f
let t0s = t0s2, as = as2, bs = bs2, sol = sol2
    Float64[as[end][0], bs[end][0], as[end][1], bs[end][1]] - sol(t0s[end])
end

# ╔═╡ b89706dc-35c7-4632-aedb-e024edc676ef
let t0s = t0s3, as = as3, bs = bs3, sol = sol3
    Float64[as[end][0], bs[end][0], as[end][1], bs[end][1]] - sol(t0s[end])
end

# ╔═╡ f75df121-e82e-426a-b69d-10c1a1dbe11f
let t0s = t0s4, as = as4, bs = bs4, sol = sol4
    Float64[as[end][0], bs[end][0], as[end][1], bs[end][1]] - sol(t0s[end])
end

# ╔═╡ 47228157-19cf-496b-98a3-6c438870763f
let t0s = t0s5, as = as5, bs = bs5, sol = sol5
    Float64[as[end][0], bs[end][0], as[end][1], bs[end][1]] - sol(t0s[end])
end

# ╔═╡ d192f1c2-22f8-45d4-bb4f-19c7a183afeb
let t0s = t0s6, as = as6, bs = bs6, sol = sol6
    Float64[as[end][0], bs[end][0], as[end][1], bs[end][1]] - sol(t0s[end])
end

# ╔═╡ 564c8e19-890a-46d1-ab46-49a51c5c7d8d
let t0s = t0s1, as = as1, bs = bs1, sol = sol1
    ts = range(t0s[1], t0s[end], 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    for n in eachindex(t0s, as, bs)
        scatter!(
            pl1,
            t0s[n:n],
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
        )
        plot!(
            pl1,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    for n in eachindex(t0s, as, bs)
        scatter!(pl2, t0s[n:n], t -> as[n](t - t0s[n]), color = :red, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> as[n](t - t0s[n]),
            color = :red,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
        scatter!(pl2, t0s[n:n], t -> bs[n](t - t0s[n]), color = :blue, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> bs[n](t - t0s[n]),
            color = :blue,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    plot(pl1, pl2)
end

# ╔═╡ 628d75e6-7af4-44a9-8f8c-1e183272bb9b
let t0s = t0s2, as = as2, bs = bs2, sol = sol2
    ts = range(t0s[1], t0s[end], 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    for n in eachindex(t0s, as, bs)
        scatter!(
            pl1,
            t0s[n:n],
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
        )
        plot!(
            pl1,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    for n in eachindex(t0s, as, bs)
        scatter!(pl2, t0s[n:n], t -> as[n](t - t0s[n]), color = :red, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> as[n](t - t0s[n]),
            color = :red,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
        scatter!(pl2, t0s[n:n], t -> bs[n](t - t0s[n]), color = :blue, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> bs[n](t - t0s[n]),
            color = :blue,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    plot(pl1, pl2)
end

# ╔═╡ 6a25736f-8219-43ee-aa7f-8bf6bfd31b74
let t0s = t0s3, as = as3, bs = bs3, sol = sol3
    ts = range(t0s[1], t0s[end], 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    for n in eachindex(t0s, as, bs)
        scatter!(
            pl1,
            t0s[n:n],
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
        )
        plot!(
            pl1,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    for n in eachindex(t0s, as, bs)
        scatter!(pl2, t0s[n:n], t -> as[n](t - t0s[n]), color = :red, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> as[n](t - t0s[n]),
            color = :red,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
        scatter!(pl2, t0s[n:n], t -> bs[n](t - t0s[n]), color = :blue, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> bs[n](t - t0s[n]),
            color = :blue,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    plot(pl1, pl2)
end

# ╔═╡ 9b3e0bf6-fc08-4f84-a694-267e5b2ad40a
let t0s = t0s4, as = as4, bs = bs4, sol = sol4
    ts = range(t0s[1], t0s[end], 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    for n in eachindex(t0s, as, bs)
        scatter!(
            pl1,
            t0s[n:n],
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
        )
        plot!(
            pl1,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    for n in eachindex(t0s, as, bs)
        scatter!(pl2, t0s[n:n], t -> as[n](t - t0s[n]), color = :red, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> as[n](t - t0s[n]),
            color = :red,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
        scatter!(pl2, t0s[n:n], t -> bs[n](t - t0s[n]), color = :blue, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> bs[n](t - t0s[n]),
            color = :blue,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    plot(pl1, pl2)
end

# ╔═╡ 7ef4a83a-4a01-4634-80e8-63c17080baab
let t0s = t0s5, as = as5, bs = bs5, sol = sol5
    ts = range(t0s[1], t0s[end], 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    for n in eachindex(t0s, as, bs)
        scatter!(
            pl1,
            t0s[n:n],
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
        )
        plot!(
            pl1,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    for n in eachindex(t0s, as, bs)
        scatter!(pl2, t0s[n:n], t -> as[n](t - t0s[n]), color = :red, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> as[n](t - t0s[n]),
            color = :red,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
        scatter!(pl2, t0s[n:n], t -> bs[n](t - t0s[n]), color = :blue, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> bs[n](t - t0s[n]),
            color = :blue,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    plot(pl1, pl2)
end

# ╔═╡ d87ace86-cb18-4148-af9a-e9c268db50c3
let t0s = t0s6, as = as6, bs = bs6, sol = sol6
    ts = range(t0s[1], t0s[end], 1000)

    pl1 = plot(
        sol,
        idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2),
        tspan = (ts[1], ts[end]),
        label = "|Q(ξ)|",
    )
    for n in eachindex(t0s, as, bs)
        scatter!(
            pl1,
            t0s[n:n],
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
        )
        plot!(
            pl1,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> abs(as[n](t - t0s[n]) + im * bs[n](t - t0s[n])),
            color = n,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2)],
        tspan = (ts[1], ts[end]),
        label = ["Re(Q)" "Im(Q)"],
        color = [:red :blue],
    )
    for n in eachindex(t0s, as, bs)
        scatter!(pl2, t0s[n:n], t -> as[n](t - t0s[n]), color = :red, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> as[n](t - t0s[n]),
            color = :red,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
        scatter!(pl2, t0s[n:n], t -> bs[n](t - t0s[n]), color = :blue, label = "")
        plot!(
            pl2,
            range(max(t0s[n] - 0.5, 0.0), t0s[n] + 0.5, 100),
            t -> bs[n](t - t0s[n]),
            color = :blue,
            label = "",
            linestyle = :dash,
            linewidth = 2,
        )
    end

    plot(pl1, pl2)
end

# ╔═╡ 2f37eb97-ffb4-4c5b-a130-b4fbcc666d21
md"""
# Parameter plots
"""

# ╔═╡ 94e9afd1-d5c5-466d-94ba-215023b8e4f3
@bind ω_slider Slider(range(1.0, 3.0, 1000), default = 1.0, show_value = true)

# ╔═╡ e760161a-370d-4907-841b-6382d8bf3122
@bind κ_slider Slider(range(0.25, 2.0, 1000), default = 0.49323, show_value = true)

# ╔═╡ c9197241-2c24-4653-a021-e2c21eee0215
@bind σ_slider Slider(range(1.0, 2.5, 1000), default = 2.3, show_value = true)

# ╔═╡ 0d173477-5b67-4a3b-8777-b0a9f6443e5b
@bind μ_slider Slider(range(0.5, 2.5, 1000), default = 0.78308, show_value = true)

# ╔═╡ 727fd02a-7018-4e29-bc3c-057a0b9c7751
parameters_slider =
    (d = 1.0, ω = ω_slider, κ = κ_slider, σ = σ_slider, μ = μ_slider, ϵ = 0.0, δ = 0.0)

# ╔═╡ 5edc2bb1-b4c0-4777-b0d5-69e013d2ca72
let parameters =
        [parameters1, parameters2, parameters3, parameters4, parameters5, parameters6]
    (
        extrema(getfield.(parameters, :κ)),
        extrema(getfield.(parameters, :ω)),
        extrema(getfield.(parameters, :σ)),
        extrema(getfield.(parameters, :μ)),
    )

end

# ╔═╡ 88237c11-e403-4d60-806e-7019945e5627
sol_slider = let p = parameters_slider
    prob = ODEProblem(f!, [p.μ, 0.0, 0.0, 0.0], tspan, p)
    solve(prob, abstol = 1e-9, reltol = 1e-9)
end;

# ╔═╡ ea65b7b8-7600-45af-bfaa-fc095124ec92
let sol = sol_slider
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])
    plot(pl1, pl2)
end

# ╔═╡ 1f63fe06-770a-47c4-8382-f1d9b4911953
md"""
# The interval $(ξ₀, ∞)$
"""

# ╔═╡ 286d2d4d-afc4-47b7-9fdc-0863593cf757
setprecision(Arb, 64)

# ╔═╡ 4d3ac31c-b4b5-4caf-ac0b-9c36563cf6ef
parameters1

# ╔═╡ 7e60fa9a-620f-459f-9a1d-8f6284587b97
function P(d, κ, ω, σ, ϵ, ξ)
    a = (1 / σ + im * ω / κ) / 2
    b = d / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    return hypgeom_u(a, b, z)
end

# ╔═╡ 0a56fa1e-37c4-4295-a206-4b71da9a06e4
function P_dξ(d, κ, ω, σ, ϵ, ξ)
    a = (1 / σ + im * ω / κ) / 2
    b = d / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    return hypgeom_u_dz(a, b, z) * (-im * κ / (1 - im * ϵ) * ξ)
end

# ╔═╡ 059e1666-2c5b-4ce9-abf7-9af61bc21e96
function E(d, κ, ω, σ, ϵ, ξ)
    a = (1 / σ + im * ω / κ) / 2
    b = d / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    return exp(z) * hypgeom_u(b - a, b, -z)
end

# ╔═╡ 8ca25816-f74c-4c42-b63b-a237dae56431
function W(d, κ, ω, σ, ϵ, ξ)
    a = (1 / σ + im * ω / κ) / 2
    b = d / 2
    z = -im * κ / (1 - im * ϵ) * ξ^2 / 2
    return -im * κ / (1 - im * ϵ) *
           exp(sign(real(z)) * im * (π * (b - a))) *
           ξ *
           z^-b *
           exp(z)
end

# ╔═╡ 4fdab699-bccb-4da0-b807-eef170d9f4c3
function K(d, κ, ω, σ, ϵ, ξ, η)
    if η <= ξ
        return -1 / (1 - im * ϵ) * P(d, κ, ω, σ, ϵ, ξ) * E(d, κ, ω, σ, ϵ, η) /
               W(d, κ, ω, σ, ϵ, η)
    elseif ξ <= η
        return -1 / (1 - im * ϵ) * E(d, κ, ω, σ, ϵ, ξ) * P(d, κ, ω, σ, ϵ, η) /
               W(d, κ, ω, σ, ϵ, η)
    else
        error("What to do in this case?")
    end
end

# ╔═╡ 5ff57f9a-7494-49a7-b663-d0e6c79e58ff
let p = parameters1
    (
        P(p.d, p.κ, p.ω, p.σ, p.ϵ, 30.0),
        E(p.d, p.κ, p.ω, p.σ, p.ϵ, 30.0),
        W(p.d, p.κ, p.ω, p.σ, p.ϵ, 30.0),
        K(p.d, p.κ, p.ω, p.σ, p.ϵ, 30.0, 35.0),
    )
end

# ╔═╡ 7ca4cbcb-a07d-4fd8-81c4-c5b3af0f67ce
let p = parameters1
    ξs = range(0.0, 30.0, 1000)
    res = P.(p.d, p.κ, p.ω, p.σ, p.ϵ, ξs)
    plot(ξs, [real.(res) imag.(res)])
end

# ╔═╡ d359870e-7747-4f12-b7f9-cb057bf2999e
md"""
# Combining $(0, ξ₀)$ and $(ξ₀, ∞)$
We plot the numerical solution on $(0, ξ₀)$ and the linear part of the solution on $(ξ₀, ∞)$.
"""

# ╔═╡ 1c891260-e705-45a4-9ec1-188a4a5ec5ad
let sol = sol1, p = parameters1
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(
        sol,
        idxs = [(0, 1), (0, 2), (0, 3), (0, 4)],
        label = ["Re(Q)" "Im(Q)" "Re(Q')" "Im(Q')"],
    )

    λ = (sol[end][1] + im * sol[end][2]) / P.(p.d, p.κ, p.ω, p.σ, p.ϵ, sol.t[end])

    ξs = range(0.1sol.t[end], 2sol.t[end], 100)
    res = λ .* P.(p.d, p.κ, p.ω, p.σ, p.ϵ, ξs)
    dres = [λ .* P.(p.d, p.κ, p.ω, p.σ, p.ϵ, ArbSeries((ξ, 1)))[1] for ξ in ξs]

    plot!(pl1, ξs, abs.(res), label = "λP")
    plot!(
        pl2,
        ξs,
        [real.(res) imag.(res) real.(dres) imag.(dres)],
        color = [1 2 3 4],
        linestyle = :dash,
        label = ["Re(Q)" "Im(Q)" "Re(Q')" "Im(Q')"],
        m = :circle,
        ms = 1,
    )

    plot(pl1, pl2)
end

# ╔═╡ addd62fc-3abb-4d3a-b53e-54c57081d4f4
let sol = sol2, p = parameters2
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])

    λ = (sol[end][1] + im * sol[end][2]) / P.(p.d, p.κ, p.ω, p.σ, p.ϵ, sol.t[end])

    ξs = range(0.1sol.t[end], 2sol.t[end], 100)
    res = λ .* P.(p.d, p.κ, p.ω, p.σ, p.ϵ, ξs)

    plot!(pl1, ξs, abs.(res), label = "λP")
    plot!(
        pl2,
        ξs,
        [real.(res) imag.(res)],
        color = [1 2],
        linestyle = :dash,
        label = ["Re(λP)" "Im(λP)"],
    )

    plot(pl1, pl2)
end

# ╔═╡ eece93a9-24b1-40d7-993c-5dfc5acfd3fc
let sol = sol3, p = parameters3
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])

    λ = (sol[end][1] + im * sol[end][2]) / P.(p.d, p.κ, p.ω, p.σ, p.ϵ, sol.t[end])

    ξs = range(0.5sol.t[end], 2sol.t[end], 100)
    res = λ .* P.(p.d, p.κ, p.ω, p.σ, p.ϵ, ξs)

    plot!(pl1, ξs, abs.(res), label = "λP")
    plot!(
        pl2,
        ξs,
        [real.(res) imag.(res)],
        color = [1 2],
        linestyle = :dash,
        label = ["Re(λP)" "Im(λP)"],
    )

    plot(pl1, pl2)
end

# ╔═╡ 33b8b346-63c3-4104-8abc-5e43b48248d4
let sol = sol4, p = parameters4
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])

    λ = (sol[end][1] + im * sol[end][2]) / P.(p.d, p.κ, p.ω, p.σ, p.ϵ, sol.t[end])

    ξs = range(0.1sol.t[end], 2sol.t[end], 100)
    res = λ .* P.(p.d, p.κ, p.ω, p.σ, p.ϵ, ξs)

    plot!(pl1, ξs, abs.(res), label = "λP")
    plot!(
        pl2,
        ξs,
        [real.(res) imag.(res)],
        color = [1 2],
        linestyle = :dash,
        label = ["Re(λP)" "Im(λP)"],
    )

    plot(pl1, pl2)
end

# ╔═╡ 9c8e04d8-5d4a-43c2-88bd-db9cf74b47c2
let sol = sol5, p = parameters5
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])

    λ = (sol[end][1] + im * sol[end][2]) / P.(p.d, p.κ, p.ω, p.σ, p.ϵ, sol.t[end])

    ξs = range(0.15sol.t[end], 2sol.t[end], 100)
    res = λ .* P.(p.d, p.κ, p.ω, p.σ, p.ϵ, ξs)

    plot!(pl1, ξs, abs.(res), label = "λP")
    plot!(
        pl2,
        ξs,
        [real.(res) imag.(res)],
        color = [1 2],
        linestyle = :dash,
        label = ["Re(λP)" "Im(λP)"],
    )

    plot(pl1, pl2)
end

# ╔═╡ 79001c21-9e55-40c8-92c8-e2cdb2ec0a74
let sol = sol6, p = parameters6
    pl1 = plot(sol, idxs = ((t, x, y) -> (t, abs(x + im * y)), 0, 1, 2), label = "|Q(ξ)|")
    pl2 = plot(sol, idxs = [(0, 1), (0, 2)], label = ["Re(Q)" "Im(Q)"])

    λ = (sol[end][1] + im * sol[end][2]) / P.(p.d, p.κ, p.ω, p.σ, p.ϵ, sol.t[end])

    ξs = range(0.15sol.t[end], 2sol.t[end], 100)
    res = λ .* P.(p.d, p.κ, p.ω, p.σ, p.ϵ, ξs)

    plot!(pl1, ξs, abs.(res), label = "λP")
    plot!(
        pl2,
        ξs,
        [real.(res) imag.(res)],
        color = [1 2],
        linestyle = :dash,
        label = ["Re(λP)" "Im(λP)"],
    )

    plot(pl1, pl2)
end

# ╔═╡ a13e74ff-1d7f-4c42-bd35-4bbfa0e94c4c
md"""
## Phrase as zero-finding problem
For this we want to set up two functions, one for computing the solution starting from zero and one for computing it using the asymptotic expansion from infinity. We then look the their difference at some point inbetween. These functions should take as argument the parameters that we want to vary to find a solution.
"""

# ╔═╡ f64b2949-a37d-4472-9067-71506bbba917
parameters1

# ╔═╡ 09759727-9fd1-4ede-87c7-1cc7a95f87e9
d = 1.0

# ╔═╡ 3c1a91e8-fcbe-485d-98f8-1dc2d192e083
ω = 1.0

# ╔═╡ a202f84e-e242-4a74-8441-18e14f166a2d
ϵ = 0.0

# ╔═╡ 3b43c1fe-9ace-4f73-8e8c-9c5daae61c3f
δ = 0.0

# ╔═╡ 5f24ac99-f57a-4272-8562-e53aa1a29262
parameters1

# ╔═╡ 74ea6a28-580a-4e13-8b99-6e8209b04151
function f_zero(κ, μ; d = 1.0, ω = 1.0, σ = 2.3, ϵ = 0.0, δ = 0.0)
    p = (; d, ω, κ, σ, μ, ϵ, δ)
    prob = ODEProblem(f!, [p.μ, 0.0, 0.0, 0.0], tspan, p)
    sol = solve(prob, abstol = 1e-9, reltol = 1e-9)
    return (sol[end][1] + im * sol[end][2], sol[end][3] + im * sol[end][4])
end

# ╔═╡ 6fb90e92-0c48-4a9d-9760-fbe5bca23c0a
function f_inf(κ, μ; d = 1.0, ω = 1.0, σ = 2.3, ϵ = 0.0, δ = 0.0)
    res = P(d, κ, ω, σ, ϵ, tspan[2])
    dres = P_dξ(d, κ, ω, σ, ϵ, tspan[2])
    return (res, dres)
end

# ╔═╡ faca5018-bf75-4d91-bc7d-fb620649b2f3
q1 = f_zero(parameters1.κ, parameters1.μ)

# ╔═╡ 128ab467-f727-4193-ac33-a6b3907c5874
q2 = f_inf(parameters1.κ, parameters1.μ)

# ╔═╡ 79633b5d-84bf-460d-9940-e3d3916e865e
γ = q1[1] / q2[1]

# ╔═╡ 1f81129f-408e-4d3c-bdf8-dbd3a939a760
γ .* q2 .- q1

# ╔═╡ f7bd05de-5254-4e51-adc3-9c330b9fa4d5
md"""
Create function we want to find zero of
"""

# ╔═╡ 004cd348-096b-4461-8df2-e7580ea687f0
function f(κ, μ)
    r11, r12 = f_zero(κ, μ)
    r21, r22 = f_inf(κ, μ)

    γ = r11 / r21
    return r12 - γ * r22
end

# ╔═╡ a9ceedfe-893c-4b87-a5bc-d9500565cd8f
f(parameters1.κ, parameters1.μ)

# ╔═╡ d95c5b9a-b765-4577-a9f2-0be449eba546
parameters1.κ, parameters1.μ

# ╔═╡ 1fb6b323-ac24-4fe8-b32a-16bba5229a92
zero_sol = nlsolve([parameters1.κ, parameters1.μ], xtol = 1e-12, ftol = 1e-12) do inp
    κ, μ = inp
    res = f(κ, μ)
    return [real(res), imag(res)]
end

# ╔═╡ 8d46b70b-feea-49e7-8d5c-fe75ce61c44c
f(zero_sol.zero[1], zero_sol.zero[2])

# ╔═╡ 32315591-1769-4ccf-b038-4dcc259bcd23
plot_box_1, plot_box_2, plot_box_3 =
    let κ = zero_sol.zero[1], μ = zero_sol.zero[2], dd = 1e-5
        κs = range(κ - dd, κ + dd, 50)
        μs = range(μ - dd, μ + dd, 50)
        res = f.(κs, μs')
        m_real = maximum(abs, real.(res))
        m_imag = maximum(abs, imag.(res))
        p1 = heatmap(
            κs,
            μs,
            real.(res),
            seriescolor = :delta,
            clims = (-m_real, m_real),
            cbar = :right,
        )
        p2 = heatmap(κs, μs, imag.(res), seriescolor = :delta, clims = (-m_imag, m_imag))
        p3 = heatmap(κs, μs, abs.(res), seriescolor = :viridis)
        p1, p2, p3
    end

# ╔═╡ 4a991485-e71a-4f5a-92a1-d367fc6aee9a
plot_box_1

# ╔═╡ 540d425a-f6fc-404f-bad1-2cdba7f9944e
plot_box_2

# ╔═╡ 08343fa9-0def-4ddd-9fc6-21260727a0ad
plot_box_3

# ╔═╡ ea3e2382-d22a-4af3-a4a3-e55ffa9cf351
md"""
# Typos in thesis
While reading through things I have found a few things that seems to be mistakes.
1. In the real and imaginary parts of the equation given right after (3.1) it should be $ϵb''$ and $ϵa''$ instead of a single derivative that is written now.
2. Figure 4 and 7 look different from the ones I have generated.
3. Last sentence of page 4: "We denote by $P$ and $E$ two linearly independent solutions to (1.4)". I believe it should be (1.7) and not (1.4)?
"""

# ╔═╡ Cell order:
# ╠═60dfdc80-ba5e-11ed-16f4-ade6c584dc7b
# ╠═374e1504-4ea4-4c24-be24-ea6a85b0a11e
# ╟─d3d1fa98-770d-4f14-bfa4-ca0a28a35563
# ╟─f57069eb-ece9-4744-aa10-31edfc25e63b
# ╠═15e38430-ca56-4b5b-90c0-0eddc6725329
# ╟─87cce4db-37c5-43f3-99b6-70467486633f
# ╠═a79b3baa-ffa8-440a-ab2d-1d7c2d068b71
# ╟─8042cafa-5e7c-49bb-b646-d0f2a92bf710
# ╠═dc647047-94e6-4670-a98e-fc07b141f0ea
# ╠═d0b548ed-ca9f-4d0a-9eb0-e3a6fdfcc689
# ╠═268e303d-b41f-4836-9a83-ae1cb0c94eec
# ╠═155768b4-5daf-457d-b0d7-c90e37dd52cd
# ╠═23dcf903-caf8-47f5-a7ff-38b5df47c6bc
# ╠═fd138fb8-b459-4888-b12d-89a838c710ff
# ╟─9c8042a5-0a86-457d-9de7-f2e059c5764c
# ╠═47b35421-7db1-4738-a20d-d79b4801cd36
# ╠═7394e042-228a-4915-abd4-dc00135b2f8a
# ╠═07763f0c-23eb-476c-a818-3767d72df5fc
# ╠═e38e45bd-cb13-4fbb-a4e0-e80aff2c48c2
# ╠═b2a0efc4-67ad-4fe7-8813-93d83bb79805
# ╠═75640143-513d-42bc-af58-55858cc4086a
# ╠═b4299bd3-6a9c-4216-83c3-34495feed058
# ╟─03f2f89c-8a8e-45d1-9dd5-05d12ad6e3a6
# ╟─27685362-25d8-4f4a-a5bf-7787d9760aaf
# ╟─85fc701b-c12e-4d0a-a404-2c214e23d979
# ╟─150d8a06-7c05-40d5-a03d-5ebff8b04558
# ╟─89c7946b-3d56-41d0-a9dc-170061f246be
# ╟─6b764bdc-cc98-4411-a021-0fea807b7a11
# ╟─ff26580d-c378-46b9-9058-7d03ae67b273
# ╟─b04575eb-25a5-4839-9ef6-213429e7485a
# ╟─e8f0a6e2-81c9-44c9-8b88-6bb44c42d1d6
# ╟─4780850f-e1d0-434f-a25a-e4501a5fef84
# ╟─1828428f-c5e8-4af5-84e2-42a504eaa386
# ╟─eebe6cd5-2a1b-484b-acf8-ea117ee6ef3c
# ╟─e78dc117-66a7-4ea0-a487-0d8591d60107
# ╟─e34216ff-716c-42f1-9bb2-f66117da6446
# ╠═15315635-2762-475f-8e5f-05361f51a483
# ╠═9e936fa9-8953-4104-becf-aa4ffc14571d
# ╠═1eb3e568-37f6-43a4-86ae-686d0ce84bc4
# ╠═13dc0e82-1a67-44a8-9060-37dc986ba518
# ╟─6a0dce19-a6ff-4e1f-a332-bf558d0e37be
# ╠═a9e9a92f-dc53-492d-bc86-9753356c84ec
# ╠═14111c0e-89c0-41cd-92bd-cf9c080a2913
# ╠═cc684d4f-0feb-4150-a446-db664c56f678
# ╠═e0a50a81-6312-4607-89f3-5f44b8b3d568
# ╟─4ed0cb0c-485b-4d2f-b028-3aff709f37d9
# ╠═8125af3e-ee25-4b84-9279-5d67f95d1772
# ╠═c09d1ce4-9012-4cea-94a4-06a255b93033
# ╟─712552b8-73e0-4245-b12c-3ac63ed13e14
# ╟─ce107d7b-2c65-4b5c-9113-d1e293c2a846
# ╠═322e24da-bdfd-4298-92a6-651043a3c03f
# ╠═3035f690-7aff-4192-b8f2-a1e78094047c
# ╠═b39e0201-6458-4f87-bbe4-f41a5e40adb0
# ╠═568061ef-fa5d-45d7-9c59-1a0592099dbd
# ╠═3da81649-a810-45ac-b520-aa7ada1fb5db
# ╠═9b18d3eb-5679-428c-98bc-41d12936541f
# ╠═453ee689-ba8c-475c-814d-2bc120e35cc7
# ╠═2ca1a635-5cf1-4aa7-a0ed-692cd87369f1
# ╟─06470b3c-df86-48e3-a69b-9a2c27921c57
# ╟─05a54cee-5829-4cca-94a8-0a408a4c9716
# ╟─5bf7172d-84ff-4da1-9e7a-e44d7fed9368
# ╟─c9a6df8d-31c4-4e86-88e9-c4c6ba426e3d
# ╟─b62d1e45-c91f-4d5e-8e0c-9001fe0481f9
# ╟─e4a7ee1f-266b-4553-92b1-6f168f63c9c2
# ╠═3bd6308d-6a33-4a5d-9f03-d50f731f1d3e
# ╠═0eb0b383-2627-42f2-8038-0ca1241a7735
# ╠═e8f7d17b-8f6f-43ba-aec7-9fb27a5af8db
# ╟─7820cb27-0fb8-44be-a4aa-dd9dbf72b999
# ╟─ea5937b1-4ce5-4152-839e-53e8588008dc
# ╟─27e53267-e7aa-4315-b73b-b0a2ad4d0ff1
# ╠═4c79b2cc-67ca-412c-acdc-76e9c659b029
# ╠═032f778d-560e-4ea9-91fa-509a5797927b
# ╠═d267c8ae-dceb-448e-b8b6-33866280cfa5
# ╠═d7a56d9a-1349-4a0b-bebf-72c9d0afe086
# ╠═077462c9-4806-4b22-8cfc-4008d72ff819
# ╠═dde47573-c285-4e10-9b17-2aae6bd8f125
# ╠═381b22aa-87dd-4805-b5b6-d16e2e0e6d8c
# ╠═4b3c9483-1483-4335-b101-5a8115157718
# ╠═ed7ad5d0-85e3-4994-a073-0c3c78c6cd7f
# ╠═b89706dc-35c7-4632-aedb-e024edc676ef
# ╠═f75df121-e82e-426a-b69d-10c1a1dbe11f
# ╠═47228157-19cf-496b-98a3-6c438870763f
# ╠═d192f1c2-22f8-45d4-bb4f-19c7a183afeb
# ╟─564c8e19-890a-46d1-ab46-49a51c5c7d8d
# ╟─628d75e6-7af4-44a9-8f8c-1e183272bb9b
# ╟─6a25736f-8219-43ee-aa7f-8bf6bfd31b74
# ╟─9b3e0bf6-fc08-4f84-a694-267e5b2ad40a
# ╟─7ef4a83a-4a01-4634-80e8-63c17080baab
# ╟─d87ace86-cb18-4148-af9a-e9c268db50c3
# ╟─2f37eb97-ffb4-4c5b-a130-b4fbcc666d21
# ╠═94e9afd1-d5c5-466d-94ba-215023b8e4f3
# ╠═e760161a-370d-4907-841b-6382d8bf3122
# ╠═c9197241-2c24-4653-a021-e2c21eee0215
# ╠═0d173477-5b67-4a3b-8777-b0a9f6443e5b
# ╟─727fd02a-7018-4e29-bc3c-057a0b9c7751
# ╟─5edc2bb1-b4c0-4777-b0d5-69e013d2ca72
# ╟─88237c11-e403-4d60-806e-7019945e5627
# ╟─ea65b7b8-7600-45af-bfaa-fc095124ec92
# ╟─1f63fe06-770a-47c4-8382-f1d9b4911953
# ╠═286d2d4d-afc4-47b7-9fdc-0863593cf757
# ╠═4d3ac31c-b4b5-4caf-ac0b-9c36563cf6ef
# ╠═7e60fa9a-620f-459f-9a1d-8f6284587b97
# ╠═0a56fa1e-37c4-4295-a206-4b71da9a06e4
# ╠═059e1666-2c5b-4ce9-abf7-9af61bc21e96
# ╠═8ca25816-f74c-4c42-b63b-a237dae56431
# ╠═4fdab699-bccb-4da0-b807-eef170d9f4c3
# ╠═5ff57f9a-7494-49a7-b663-d0e6c79e58ff
# ╠═7ca4cbcb-a07d-4fd8-81c4-c5b3af0f67ce
# ╟─d359870e-7747-4f12-b7f9-cb057bf2999e
# ╟─1c891260-e705-45a4-9ec1-188a4a5ec5ad
# ╟─addd62fc-3abb-4d3a-b53e-54c57081d4f4
# ╟─eece93a9-24b1-40d7-993c-5dfc5acfd3fc
# ╟─33b8b346-63c3-4104-8abc-5e43b48248d4
# ╟─9c8e04d8-5d4a-43c2-88bd-db9cf74b47c2
# ╟─79001c21-9e55-40c8-92c8-e2cdb2ec0a74
# ╟─a13e74ff-1d7f-4c42-bd35-4bbfa0e94c4c
# ╠═f64b2949-a37d-4472-9067-71506bbba917
# ╠═09759727-9fd1-4ede-87c7-1cc7a95f87e9
# ╠═3c1a91e8-fcbe-485d-98f8-1dc2d192e083
# ╠═a202f84e-e242-4a74-8441-18e14f166a2d
# ╠═3b43c1fe-9ace-4f73-8e8c-9c5daae61c3f
# ╠═5f24ac99-f57a-4272-8562-e53aa1a29262
# ╠═74ea6a28-580a-4e13-8b99-6e8209b04151
# ╠═6fb90e92-0c48-4a9d-9760-fbe5bca23c0a
# ╠═faca5018-bf75-4d91-bc7d-fb620649b2f3
# ╠═128ab467-f727-4193-ac33-a6b3907c5874
# ╠═79633b5d-84bf-460d-9940-e3d3916e865e
# ╠═1f81129f-408e-4d3c-bdf8-dbd3a939a760
# ╟─f7bd05de-5254-4e51-adc3-9c330b9fa4d5
# ╠═004cd348-096b-4461-8df2-e7580ea687f0
# ╠═a9ceedfe-893c-4b87-a5bc-d9500565cd8f
# ╠═d95c5b9a-b765-4577-a9f2-0be449eba546
# ╠═1fb6b323-ac24-4fe8-b32a-16bba5229a92
# ╠═8d46b70b-feea-49e7-8d5c-fe75ce61c44c
# ╠═32315591-1769-4ccf-b038-4dcc259bcd23
# ╠═4a991485-e71a-4f5a-92a1-d367fc6aee9a
# ╠═540d425a-f6fc-404f-bad1-2cdba7f9944e
# ╠═08343fa9-0def-4ddd-9fc6-21260727a0ad
# ╟─ea3e2382-d22a-4af3-a4a3-e55ffa9cf351
