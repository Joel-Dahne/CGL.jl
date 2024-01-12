"""
    CGLBranch

This is a self-contained submodule for computing the branches
corresponding to Figures 3.1 and 3.6 in
https://doi.org/10.1002/cpa.3006.

A figure similar to Figure 3.1 can be produced with
```
λ = CGL.CGLBranch.Params()
μκs = CGL.CGLBranch.sverak_initial.(1:8, (λ,))
brs = Folds.map(μκs) do (μ, κ)
    CGL.CGLBranch.branch(μ, κ, λ)
end
pl = plot()
foreach(br -> plot!(pl, br), brs)
pl
```

For Figure 3.6 you would get
```
λ = CGL.CGLBranch.Params(0.0, 3, 1.0, 1.0, 0.0, 30.0)
μκs = CGL.CGLBranch.sverak_initial.(1:5, (λ,))
brs = Folds.map(μκs) do (μ, κ)
    CGL.CGLBranch.branch(μ, κ, λ)
end
pl = plot()
foreach(br -> plot!(pl, br), brs)
pl
```
However this doesn't work nearly as well and only manages to fully
track one of the branches.
"""
module CGLBranch

using BifurcationKit
using DifferentialEquations
using NLsolve
using StaticArrays

export Params, sverak_initial, branch_F, branch

@kwdef struct Params
    ϵ::Float64 = 0.0
    d::Int = 1
    ω::Float64 = 1.0
    σ::Float64 = 2.3
    δ::Float64 = 0.0
    ξ₁::Float64 = 30.0
end

rising(x, n::Integer) = prod(i -> x + i, 0:n-1, init = one(x))

function U(a, b, z)
    N = 20
    res = sum(0:N-1) do k
        rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-z)^k)
    end
    return z^-a * res
end

U_dz(a, b, z) = -a * U(a + 1, b + 1, z)

function abc(κ, λ)
    (; d, ω, σ, ϵ) = λ
    a = (1 / σ + im * ω / κ) / 2
    b = oftype(a, d) / 2
    c = -im * κ / (1 - im * ϵ) / 2
    return a, b, c
end

function P(ξ, κ, λ)
    a, b, c = abc(κ, λ)

    return U(a, b, c * ξ^2)
end

function P_dξ(ξ, κ, λ)
    a, b, c = abc(κ, λ)
    z_dξ = 2c * ξ

    return U_dz(a, b, c * ξ^2) * z_dξ
end

function E(ξ, κ, λ)
    a, b, c = abc(κ, λ)

    z = c * ξ^2

    return exp(z) * U(b - a, b, -z)
end

function E_dξ(ξ, κ, λ)
    a, b, c = abc(κ, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ

    return exp(z) * (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ
end

function W(ξ, κ, λ)
    (; ϵ) = λ

    a, b, c = abc(κ, λ)

    z = c * ξ^2

    return -im * κ / (1 - im * ϵ) *
           exp(sign(imag(c)) * im * (b - a) * π) *
           ξ *
           z^-b *
           exp(z)
end

function B_W(κ, λ)
    (; δ) = λ

    a, b, c = abc(κ, λ)

    return -(1 + im * δ) / (im * κ) * exp(-sign(imag(c)) * im * (b - a) * π) * c^b
end

function J_P(ξ, κ, λ)
    (; ϵ, δ) = λ

    return (1 + im * δ) / (1 - im * ϵ) * P(ξ, κ, λ) / W(ξ, κ, λ)
end

function J_E(ξ, κ, λ)
    (; ϵ, δ) = λ

    return (1 + im * δ) / (1 - im * ϵ) * E(ξ, κ, λ) / W(ξ, κ, λ)
end

function system_d1(u, (κ, λ), ξ)
    (; d, ω, σ, ϵ, δ) = λ
    a, b, α, β = u

    @fastmath begin
        a2b2σ = (a^2 + b^2)^σ

        F1 = κ * ξ * β + κ / σ * b + ω * a - a2b2σ * a + δ * a2b2σ * b
        F2 = -κ * ξ * α - κ / σ * a + ω * b - a2b2σ * b - δ * a2b2σ * a

        return SVector(α, β, (F1 - ϵ * F2) / (1 + ϵ^2), (ϵ * F1 + F2) / (1 + ϵ^2))
    end
end

function system(u, (κ, λ), ξ)
    (; d, ω, σ, ϵ, δ) = λ
    a, b, α, β = u

    @fastmath begin
        a2b2σ = (a^2 + b^2)^σ

        F1 = κ * ξ * β + κ / σ * b + ω * a - a2b2σ * a + δ * a2b2σ * b
        F2 = -κ * ξ * α - κ / σ * a + ω * b - a2b2σ * b - δ * a2b2σ * a

        if !iszero(ξ)
            F1 -= (d - 1) / ξ * (α + ϵ * β)
            F2 -= (d - 1) / ξ * (β - ϵ * α)
        end

        return SVector(α, β, (F1 - ϵ * F2) / (1 + ϵ^2), (ϵ * F1 + F2) / (1 + ϵ^2))
    end
end

function branch_F(x, λ::Params)
    μ, κ = x
    (; d, σ, ξ₁) = λ

    if λ.d == 1
        prob = ODEProblem{false}(system_d1, SVector(μ, 0, 0, 0), (zero(ξ₁), ξ₁), (κ, λ))
        sol = solve(
            prob,
            AutoVern7(Rodas5P()),
            abstol = 1e-9,
            reltol = 1e-9,
            maxiters = 4000,
            verbose = false,
        )
    else
        prob = ODEProblem{false}(system, SVector(μ, 0, 0, 0), (zero(ξ₁), ξ₁), (κ, λ))
        sol = solve(
            prob,
            AutoVern7(Rodas5P()),
            abstol = 1e-9,
            reltol = 1e-9,
            maxiters = 4000,
            verbose = false,
        )
    end

    a, b, α, β = sol[end]::SVector{4,typeof(μ)}

    order = 2 # Order of approximation to use

    if sol.retcode != ReturnCode.Success
        # IMPROVE: We don't really want to compute anything in this
        # case. But it is unclear exactly what we should return then.
        # For now we at least only use the first order approximation.
        order = 1
    end

    Q_0, dQ_0 = complex(a, b), complex(α, β)

    if order == 1 # First order approximation
        γ = Q_0 / P(ξ₁, κ, λ)

        dQ_inf = γ * P_dξ(ξ₁, κ, λ)
    elseif order == 2 # Second order approximation
        p = P(ξ₁, κ, λ)
        p_dξ = P_dξ(ξ₁, κ, λ)
        e = E(ξ₁, κ, λ)
        e_dξ = E_dξ(ξ₁, κ, λ)

        _, _, c = abc(κ, λ)

        γ = nlsolve([Q_0 / p], ftol = 1e-12) do γ
            Q_0 - (
                γ[1] * p +
                e *
                B_W(κ, λ) *
                exp(-c * ξ₁^2) *
                p *
                ξ₁^(d - 2) *
                abs(γ[1] * p)^2σ *
                γ[1] *
                p / 2c
            )
        end.zero[1]

        # Compute derivative for solution at infinity with given γ
        I_E = zero(γ)
        I_P = B_W(κ, λ) * exp(-c * ξ₁^2) * p * ξ₁^(d - 2) * abs(γ * p)^2σ * γ * p / 2c

        Q_inf = γ * p + p * I_E + e * I_P

        I_E_dξ = J_E(ξ₁, κ, λ) * abs(Q_inf)^2σ * Q_inf
        I_P_dξ = -J_P(ξ₁, κ, λ) * abs(Q_inf)^2σ * Q_inf

        dQ_inf = γ * p_dξ + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ
    else
        error("unsupported order $order")
    end

    res = dQ_0 - dQ_inf

    return [real(res), imag(res)]
end

function sverak_initial(i, λ::Params = Params())
    if λ.d == 1 && λ.σ == 2.3
        μs = [1.23204, 0.78308, 1.12389, 0.88393, 1.07969, 0.92761, 1.05707, 0.94914]
        κs = [0.85310, 0.49322, 0.34680, 0.26678, 0.21621, 0.18192, 0.15667, 0.13749]
    elseif λ.d == 3 && λ.σ == 1
        μs = [
            1.885903265965844,
            0.8398821872894153,
            1.1174528173836733,
            0.9763417524521208,
            1.0400247309702308,
        ]
        κs = [
            0.9174206192134661,
            0.32125552601396035,
            0.22690292060311984,
            0.16994395856681485,
            0.13917535359820785,
        ]
    end

    return μs[i], κs[i]
end

function branch(μ, κ, λ::Params; plot = false)
    # IMPROVE: Look at using ShootingProblem
    prob = BifurcationProblem(
        branch_F,
        [μ, κ],
        λ,
        (@lens _.ϵ),
        record_from_solution = (x, λ) -> (κ = x[2], μ = x[1]),
    )

    if λ.d == 1 && λ.σ == 2.3
        opts = ContinuationPar(
            p_min = 0.0,
            p_max = 0.07,
            dsmin = 0.00005,
            ds = 0.0001,
            dsmax = 0.0005,
            max_steps = 1500,
            detect_bifurcation = 0,
            newton_options = NewtonPar(tol = 1e-6, max_iterations = 10),
        )

        finalise_solution =
            (z, tau, step, contResult; kwargs...) -> !(tau.p < 0 && z.p < 0.02)
    elseif λ.d == 3 && λ.σ == 1
        # IMPROVE: This is not able to capture the full branches. It
        # needs more tuning or other changes.
        opts = ContinuationPar(
            p_min = 0.0,
            p_max = 0.3,
            dsmin = 0.000005,
            ds = 0.0001,
            dsmax = 0.0005,
            max_steps = 2000,
            detect_bifurcation = 0,
            newton_options = NewtonPar(tol = 1e-6, max_iterations = 10),
        )

        finalise_solution =
            (z, tau, step, contResult; kwargs...) -> !(tau.p < 0 && z.p < 0.1)
    end

    br = continuation(prob, PALC(), opts; finalise_solution, plot)

    return br
end

end # module CGLBranch
