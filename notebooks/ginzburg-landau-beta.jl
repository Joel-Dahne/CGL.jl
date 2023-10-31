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

# ╔═╡ 97f27c00-ced6-11ed-3a56-83fda833255c
begin
    using Pkg, Revise
    Pkg.activate("..", io = devnull)
    using Arblib
    using DifferentialEquations
    using Folds
    using CGL
    using NLsolve
    using Plots
    using PlutoUI
    using StaticArrays
    import IntervalArithmetic
end

# ╔═╡ 69c032a8-f220-49cd-a388-2551fdd5c045
setprecision(Arb, 128)

# ╔═╡ c159ad00-9986-4a6b-b273-b64f0ace2391
params = gl_params(1, 1.0, 2.3, 0.0, 0.0)

# ╔═╡ 57e417ce-6feb-423d-858a-9c06121b265e
tspan = (0.0, 10.0)

# ╔═╡ 2db36a21-6a3d-462f-b3bd-b5829dda096b
κ_start, μ_start = 0.49323, 0.78308

# ╔═╡ 67e866d5-b54e-4b14-9b4c-28746d5dfedf
function f_zero(κ, μ; return_sol = false)
    prob = ODEProblem(gl_equation_real, SVector(μ, 0.0, 0.0, 0.0), tspan, (params, κ))
    sol = solve(prob, abstol = 1e-9, reltol = 1e-9)
    return_sol && return sol
    return SVector(sol[end][1] + im * sol[end][2], sol[end][3] + im * sol[end][4])
end

# ╔═╡ 9f312bde-9d85-4c7f-9922-5467dc66a169
function f_inf(κ, μ)
    res = P(tspan[2], (params, κ))
    dres = P_dξ(tspan[2], (params, κ))
    return SVector(res, dres)
end

# ╔═╡ b4fec6a0-3594-41d0-b79b-aa926adc4a30
function f(κ, μ)
    r11, r12 = f_zero(κ, μ)
    r21, r22 = f_inf(κ, μ)

    γ = r11 / r21
    return r12 - γ * r22
end

# ╔═╡ c6ca0713-489f-411a-a76a-15222507ce0e
sol_params = nlsolve([κ_start, μ_start], xtol = 1e-12, ftol = 1e-12) do inp
    κ, μ = inp
    res = f(κ, μ)
    return [real(res), imag(res)]
end

# ╔═╡ 5398cb24-22a1-4c9e-8902-5e4f3a0981b3
κ, μ = sol_params.zero

# ╔═╡ 5446937a-ee72-4492-88a3-320010e3f514
f(κ, μ)

# ╔═╡ 0a23fc32-3bae-4f25-a170-440df1af8405
κs, μs, res = let δκ = 1e-5, δμ = 1e-5, N = 50
    κs = range(κ - δκ, κ + δκ, N)
    μs = range(μ - δμ, μ + δμ, N)
    res = f.(κs, μs')

    κs, μs, res
end

# ╔═╡ a7a0c16c-5520-4ade-a229-cdd8f6785fc8
plot_box_1, plot_box_2, plot_box_3 = let
    # These are used to make the coloring symmetric around zero	
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

# ╔═╡ 4ab79013-e35b-49df-bf51-ba12d6061ddc
plot_box_1

# ╔═╡ d77a49f4-330a-4f9c-9c27-a8a837e28558
plot_box_2

# ╔═╡ 88834e1d-98ba-4ac6-ab0c-b9301a93159e
plot_box_3

# ╔═╡ a214b50e-3ae6-4ebe-ae65-9d6931357e89
κ_face_1 = imag.(f.(κs[begin], μs))

# ╔═╡ 3676dabd-83f1-4332-9afa-d241a7877c60
κ_face_2 = imag.(f.(κs[end], μs))

# ╔═╡ 4abbac2a-0dca-4022-820d-a0c0ff4bf50b
μ_face_1 = real.(f.(κs, μs[begin]))

# ╔═╡ 0bba527c-da1c-4364-b8d3-0d855dd88860
μ_face_2 = real.(f.(κs, μs[end]))

# ╔═╡ d803408f-9b19-4496-be99-6cfc60c9c2d6
plot(μs, [κ_face_1 κ_face_2], xlabel = "μ", labels = ["κ = κs[begin]" "κ = κs[end]"])

# ╔═╡ 48430aea-066a-46dd-ab99-d6fa77c4d1da
plot(κs, [μ_face_1 μ_face_2], xlabel = "κ", labels = ["μ = μs[begin]" "μ = μs[end]"])

# ╔═╡ 1b9dcb80-9d20-47b1-9c24-ae1fe4f50ccf
md"""
## Taylor expansion
"""

# ╔═╡ 3df308f6-7467-4a8e-b598-34782626c811
params_arb = gl_params(Arb, params)

# ╔═╡ 9afe7c7d-a013-43bc-b3af-bc07301c133d
a_series, b_series = gl_taylor_expansion_real(
    NTuple{2,Arb}[(μ, 0), (0, 0)],
    Arb(0),
    (params_arb, Arb(κ)),
    degree = 10,
)

# ╔═╡ 80c2f3ae-43f0-40f7-80d4-88a7130e076d
let ts = range(0, 1, 1000)

    pl1 = plot(ts, abs.(a_series.(ts) + im * b_series.(ts)), linestyle = :dash)
    pl2 = plot(
        ts,
        [a_series.(ts) b_series.(ts)],
        label = "",
        color = [:red :blue],
        linestyle = :dash,
    )

    plot(pl1, pl2)
end

# ╔═╡ 2a5648fc-60b0-4283-b7d5-53b41ccfba6c
prob_series = ODESeriesSecondOrderProblem(
    gl_taylor_expansion_real,
    NTuple{2,Arb}[(μ, 0), (0, 0)],
    Arb.(tspan),
    (params_arb, Arb(κ)),
)

# ╔═╡ 521cf6c2-9b65-4023-9c53-8d66fa60ffc3
sol_series = ode_series_solver(prob_series)

# ╔═╡ 670fe2cf-499c-4404-a2f7-0409a304816b
sol = f_zero(κ, μ, return_sol = true);

# ╔═╡ e47bcefe-7e3e-436b-81ab-f5198fbcb5b7
let
    pl = plot()

    plot!(pl, sol, idxs = [(0, 1), (0, 2)])

    y = [getindex.(getindex.(sol_series.u, 1), 1) getindex.(getindex.(sol_series.u, 2), 1)]
    plot!(pl, sol_series.t, y, ribbon = radius.(Arb, y), m = :circle, ms = 1)

    pl
end

# ╔═╡ 7b3af79b-0690-4ba3-bc89-94bc5ff903ca
plot(
    sol_series.t[2:end],
    radius.(Arb, getindex.(getindex.(sol_series.u, 1), 1))[2:end],
    yaxis = :log10,
)

# ╔═╡ b03da3e0-3729-4b56-b1e8-3958c8b45f82
@bind solver_degree Slider(1:40, default = 5, show_value = true)

# ╔═╡ 96928d62-ed8c-4dd4-bad3-e523fffd9b8c
@bind solver_Δt Slider(range(0, 0.5, 50)[2:end], default = 0.1, show_value = true)

# ╔═╡ b84eda81-d129-40fb-a15b-bb3621116cb2
@bind solver_prec Slider(32:8:512, default = 128, show_value = true)

# ╔═╡ fc7a97ca-b22d-45c8-bccc-b2b419ec2b21
sol_series_vary = setprecision(Arb, solver_prec) do
    prob_series_vary = ODESeriesSecondOrderProblem(
        gl_taylor_expansion_real,
        NTuple{2,Arb}[(setball(Arb, μ, Mag(0)), 0), (0, 0)],
        Arb.(tspan),
        (params_arb, setball(Arb, κ, Mag(1e-40))),
    )

    ode_series_solver(prob_series_vary; degree = solver_degree, Δt = Arb(solver_Δt))
end

# ╔═╡ afad916c-547e-4797-9de6-403b6525fc66
sol_series_vary_auton = setprecision(Arb, solver_prec) do
    prob_series_vary = ODESeriesAutonomusProblem(
        gl_taylor_expansion_real_autonomus,
        Arb[0, μ, 0, 0, 0],
        Arb.(tspan),
        (params_arb, setball(Arb, κ, Mag(1e-40))),
    )
    ode_series_solver(prob_series_vary; degree = solver_degree, Δt = Arb(solver_Δt))
end

# ╔═╡ 5423ecea-c01b-4e52-9780-66c9221268ac
let
    pl = plot(yaxis = :log10)

    plot!(
        pl,
        sol_series.t[2:end],
        radius.(Arb, getindex.(getindex.(sol_series.u, 1), 1))[2:end],
    )

    xl, yl = xlims(pl), ylims(pl)

    plot!(
        pl,
        sol_series_vary.t[2:end],
        radius.(Arb, getindex.(getindex.(sol_series_vary.u, 1), 1))[2:end],
    )

    plot!(
        pl,
        sol_series_vary_auton.t[2:end],
        radius.(Arb, getindex.(sol_series_vary_auton.u, 2))[2:end],
    )

    #plot!(pl, xlims = xl, ylims = yl)

    pl
end

# ╔═╡ fcf7349d-418b-4a2f-b559-7aa907e9a351
let
    pl = plot()

    diff_1 = getindex.(getindex.(sol_series.u, 1), 1) .- getindex.(sol.(sol_series.t), 1)
    diff_2 =
        getindex.(getindex.(sol_series_vary.u, 1), 1) .-
        getindex.(sol.(sol_series_vary.t), 1)
    diff_3 =
        getindex.(sol_series_vary_auton.u, 2) .- getindex.(sol.(sol_series_vary_auton.t), 1)

    plot!(
        pl,
        sol_series.t,
        diff_1,
        #ribbon = radius.(Arb, diff_1),
        m = :circle,
        ms = 1,
    )

    xl, yl = xlims(pl), ylims(pl)

    plot!(
        pl,
        sol_series_vary.t,
        diff_2,
        #ribbon = radius.(Arb, diff_2),
        m = :circle,
        ms = 1,
    )
    plot!(
        pl,
        sol_series_vary_auton.t,
        diff_3,
        #ribbon = radius.(Arb, diff_2),
        m = :circle,
        ms = 1,
    )

    plot!(pl, xlims = xl, ylims = (-1e-5, 1e-5))
    pl
end

# ╔═╡ b08f13e1-6788-424f-9f21-439871692eb6
@bind solver_degree_2 Slider(1:40, default = 5, show_value = true)

# ╔═╡ 1f1b18f4-d656-4aa9-9fa0-963de6d443d6
@bind solver_Δt_2 Slider(range(0, 0.5, 50)[2:end], default = 0.1, show_value = true)

# ╔═╡ a856ab98-4f39-4eec-adf6-b61a53092c52
@bind solver_prec_2 Slider(32:8:512, default = 128, show_value = true)

# ╔═╡ 6e81a5e6-2229-4eab-9667-94422b6ed521
sol_series_vary_auton2 = setprecision(Arb, solver_prec_2) do
    prob_series_vary = ODESeriesAutonomusProblem(
        gl_taylor_expansion_real_autonomus,
        Arb[0, μ, 0, 0, 0],
        Arb.(tspan),
        (gl_params(Arb, params), setball(Arb, κ, Mag(1e-30))),
    )
    ode_series_solver(
        prob_series_vary;
        degree = solver_degree_2,
        Δt = Arb(solver_Δt_2),
        skip_remainder = true,
    )
end

# ╔═╡ ed0ab333-e966-449a-9633-2c5cd0bef791
sol_series_vary_auton2_simple = setprecision(Arb, solver_prec_2) do
    prob_series_vary = ODESeriesAutonomusProblem(
        gl_taylor_expansion_real_autonomus_simple,
        Arb[0, μ, 0, 0, 0],
        Arb.(tspan),
        (gl_params(Arb, params), setball(Arb, κ, Mag(1e-30))),
    )
    ode_series_solver(
        prob_series_vary;
        degree = solver_degree_2,
        Δt = Arb(solver_Δt_2),
        skip_remainder = true,
    )
end

# ╔═╡ 9211ef97-7c73-4989-a216-b2d4129441c8
let
    pl = plot(yaxis = :log10)

    plot!(
        pl,
        sol_series_vary_auton2.t[2:end],
        radius.(Arb, sol_series_vary_auton2[:, 2][2:end]),
    )

    plot!(
        pl,
        sol_series_vary_auton2_simple.t[2:end],
        radius.(Arb, sol_series_vary_auton2_simple[:, 2][2:end]),
    )

    pl
end

# ╔═╡ d1ac4ae9-9b65-4c3c-9802-bbb4fa9da3bc
md"""
# DifferentialEquations with Arb
"""

# ╔═╡ 12c2dee8-36b5-4160-b3d0-453d130b961e
function f_zero_test(κ, u0, tspan, params, tol = 1e-9, solver = nothing)
    f = (u, p, ξ) -> begin
        gl_equation_real(u, p, ξ)
    end

    prob = ODEProblem(f, u0, tspan, (params, κ))
    sol = if isnothing(solver)
        solve(prob, abstol = tol, reltol = tol, dt = oftype(tspan[1], 0.1))
    else
        solve(prob, solver, abstol = tol, reltol = tol, dt = oftype(tspan[1], 0.1))
    end
    return sol
end

# ╔═╡ 1e130723-9551-4001-90ea-6ffc20289a9c
sol_float64 = f_zero_test(κ, sol(1), (1.0, 3.0), params);

# ╔═╡ 119540e7-53c6-4077-b834-c26678f56bb6
sol_float64.alg

# ╔═╡ b60b0d21-30aa-44b0-a063-2c9e17588fbf
Base.nextfloat(x::Arb) = begin
    x + eps(x)
end

# ╔═╡ 8a5fbfa7-09ba-4d27-9c21-1fac238b6ec6
Base.log10(x::Arb) = log(x) / log(Arb(10))

# ╔═╡ 11fd7228-8054-4a8c-aee0-92e2ff3bd95a
sol_arb = f_zero_test(Arb(κ), Arb.(sol(1)), Arb.((1.0, 3.0)), gl_params(Arb, params));

# ╔═╡ 184c7456-e118-4245-807d-b7e1172ab5d1
sol_arb_2 = f_zero_test(
    Arb(κ),
    Arb.(sol(1)),
    Arb.((1.0, 3.0)),
    gl_params(Arb, params),
    1e-4,
    AutoTsit5(Rosenbrock23()),
);

# ╔═╡ d6afc8f0-957c-44bb-9cf6-e5433a69357a
sol_arb.alg

# ╔═╡ bbb01d95-d2af-4696-b710-a864375d959f
sol_arb_2.alg

# ╔═╡ ca3c9798-75ea-4255-8de7-153e72a30c84
sol(1)

# ╔═╡ 9d510800-2a8f-4119-ba45-555f89cfa89c
sol_arb_series = setprecision(Arb, precision(Arb)) do
    prob_series = ODESeriesAutonomusProblem(
        gl_taylor_expansion_real_autonomus,
        [1.0; Arb.(sol(1))],
        Arb.((1.0, 3.0)),
        (gl_params(Arb, params), Arb(κ)),
    )
    ode_series_solver(prob_series; degree = 40, Δt = Arb(0.1), skip_remainder = true)
end

# ╔═╡ f4b10b86-2f10-46b4-b28a-1a995cb148b1
Base.Float32(x::IntervalArithmetic.Interval) = Float32(IntervalArithmetic.mid(x))

# ╔═╡ ff6bdb52-72b4-463f-9f52-a4088b88a34b
Base.Float64(x::IntervalArithmetic.Interval) = Float64(IntervalArithmetic.mid(x))

# ╔═╡ 100bf3fc-1ed7-4ba8-b0af-d13e4c779dad
Base.Float32(x::Arb) = Float32(Float64(x))

# ╔═╡ 7751f625-7362-4326-ae19-4586a1e1500c
sol_interval = f_zero_test(
    IntervalArithmetic.interval(κ),
    IntervalArithmetic.interval.(sol(1)),
    IntervalArithmetic.interval.((1.0, 3.0)),
    gl_params(IntervalArithmetic.Interval{Float64}, params),
)

# ╔═╡ d5f1c9d0-fe7b-47f4-877b-38d7af504386
let
    pl = plot()

    plot!(pl, sol_float64, idxs = [(0, 1)], linewidth = 10)

    plot!(pl, sol_arb, idxs = [(0, 1)], linewidth = 4)

    plot!(pl, sol_arb_2, idxs = [(0, 1)], linewidth = 4)

    plot!(pl, sol_arb_series.t, sol_arb_series[:, 2])

    plot!(
        pl,
        IntervalArithmetic.mid.(sol_interval.t),
        IntervalArithmetic.mid.(getindex.(sol_interval.u, 1)),
        m = :circle,
        label = "Interval",
    )

    pl
end

# ╔═╡ 6fd18525-0a90-4772-838e-6a0ca0442fbb
let
    pl = plot(yaxis = :log10)

    plot!(pl, sol_arb.t[2:end], radius.(Arb, getindex.(sol_arb.u, 1))[2:end], m = :circle)

    plot!(
        pl,
        sol_arb_2.t[2:end],
        radius.(Arb, getindex.(sol_arb_2.u, 1))[2:end],
        m = :circle,
    )

    plot!(
        pl,
        sol_arb_series.t[2:end],
        radius.(Arb, sol_arb_series[:, 2])[2:end],
        m = :square,
    )

    plot!(
        pl,
        sol_interval.t[2:end],
        IntervalArithmetic.radius.(getindex.(sol_interval.u, 1))[2:end],
        m = :diamond,
    )
end

# ╔═╡ Cell order:
# ╠═97f27c00-ced6-11ed-3a56-83fda833255c
# ╠═69c032a8-f220-49cd-a388-2551fdd5c045
# ╠═c159ad00-9986-4a6b-b273-b64f0ace2391
# ╠═57e417ce-6feb-423d-858a-9c06121b265e
# ╠═2db36a21-6a3d-462f-b3bd-b5829dda096b
# ╠═67e866d5-b54e-4b14-9b4c-28746d5dfedf
# ╠═9f312bde-9d85-4c7f-9922-5467dc66a169
# ╠═b4fec6a0-3594-41d0-b79b-aa926adc4a30
# ╠═c6ca0713-489f-411a-a76a-15222507ce0e
# ╠═5398cb24-22a1-4c9e-8902-5e4f3a0981b3
# ╠═5446937a-ee72-4492-88a3-320010e3f514
# ╠═0a23fc32-3bae-4f25-a170-440df1af8405
# ╠═a7a0c16c-5520-4ade-a229-cdd8f6785fc8
# ╠═4ab79013-e35b-49df-bf51-ba12d6061ddc
# ╠═d77a49f4-330a-4f9c-9c27-a8a837e28558
# ╠═88834e1d-98ba-4ac6-ab0c-b9301a93159e
# ╠═a214b50e-3ae6-4ebe-ae65-9d6931357e89
# ╠═3676dabd-83f1-4332-9afa-d241a7877c60
# ╠═4abbac2a-0dca-4022-820d-a0c0ff4bf50b
# ╠═0bba527c-da1c-4364-b8d3-0d855dd88860
# ╠═d803408f-9b19-4496-be99-6cfc60c9c2d6
# ╠═48430aea-066a-46dd-ab99-d6fa77c4d1da
# ╟─1b9dcb80-9d20-47b1-9c24-ae1fe4f50ccf
# ╠═3df308f6-7467-4a8e-b598-34782626c811
# ╠═9afe7c7d-a013-43bc-b3af-bc07301c133d
# ╠═80c2f3ae-43f0-40f7-80d4-88a7130e076d
# ╠═2a5648fc-60b0-4283-b7d5-53b41ccfba6c
# ╠═521cf6c2-9b65-4023-9c53-8d66fa60ffc3
# ╠═670fe2cf-499c-4404-a2f7-0409a304816b
# ╟─e47bcefe-7e3e-436b-81ab-f5198fbcb5b7
# ╟─7b3af79b-0690-4ba3-bc89-94bc5ff903ca
# ╠═b03da3e0-3729-4b56-b1e8-3958c8b45f82
# ╠═96928d62-ed8c-4dd4-bad3-e523fffd9b8c
# ╠═b84eda81-d129-40fb-a15b-bb3621116cb2
# ╠═fc7a97ca-b22d-45c8-bccc-b2b419ec2b21
# ╠═afad916c-547e-4797-9de6-403b6525fc66
# ╟─5423ecea-c01b-4e52-9780-66c9221268ac
# ╟─fcf7349d-418b-4a2f-b559-7aa907e9a351
# ╠═b08f13e1-6788-424f-9f21-439871692eb6
# ╠═1f1b18f4-d656-4aa9-9fa0-963de6d443d6
# ╠═a856ab98-4f39-4eec-adf6-b61a53092c52
# ╠═6e81a5e6-2229-4eab-9667-94422b6ed521
# ╟─ed0ab333-e966-449a-9633-2c5cd0bef791
# ╠═9211ef97-7c73-4989-a216-b2d4129441c8
# ╟─d1ac4ae9-9b65-4c3c-9802-bbb4fa9da3bc
# ╠═12c2dee8-36b5-4160-b3d0-453d130b961e
# ╠═1e130723-9551-4001-90ea-6ffc20289a9c
# ╠═119540e7-53c6-4077-b834-c26678f56bb6
# ╠═b60b0d21-30aa-44b0-a063-2c9e17588fbf
# ╠═100bf3fc-1ed7-4ba8-b0af-d13e4c779dad
# ╠═8a5fbfa7-09ba-4d27-9c21-1fac238b6ec6
# ╠═11fd7228-8054-4a8c-aee0-92e2ff3bd95a
# ╠═184c7456-e118-4245-807d-b7e1172ab5d1
# ╠═d6afc8f0-957c-44bb-9cf6-e5433a69357a
# ╠═bbb01d95-d2af-4696-b710-a864375d959f
# ╠═ca3c9798-75ea-4255-8de7-153e72a30c84
# ╠═9d510800-2a8f-4119-ba45-555f89cfa89c
# ╠═f4b10b86-2f10-46b4-b28a-1a995cb148b1
# ╠═ff6bdb52-72b4-463f-9f52-a4088b88a34b
# ╠═7751f625-7362-4326-ae19-4586a1e1500c
# ╠═d5f1c9d0-fe7b-47f4-877b-38d7af504386
# ╠═6fd18525-0a90-4772-838e-6a0ca0442fbb
