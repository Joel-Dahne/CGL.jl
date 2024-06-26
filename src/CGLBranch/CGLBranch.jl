"""
    CGLBranch

This is a self-contained submodule for computing the branches
corresponding to Figures 3.1 and 3.6 in
https://doi.org/10.1002/cpa.3006.

A figure similar to Figure 3.1 can be produced with
```
brs = tmap(1:8) do j
        CGL.CGLBranch.branch_epsilon(CGL.CGLBranch.sverak_initial(j, 1)...)
end
pl = plot()
foreach(br -> plot!(pl, br), brs)
pl
```

For Figure 3.6 you would get
```
brs3 = tmap(1:5) do j
    CGL.CGLBranch.branch_epsilon(CGL.CGLBranch.sverak_initial(j, 3)...)
end
pl = plot()
foreach(br -> plot!(pl, br), brs)
pl
```

The above examples do the continuation in `ϵ`. It is also possible to
do the continuation in `κ` instead by replacing `_epsilon` with
`_kappa` in the method.
"""
module CGLBranch

using BifurcationKit
using DifferentialEquations
using NLsolve
using StaticArrays

@kwdef struct Params
    d::Int = 1
    σ::Float64 = 2.3
    δ::Float64 = 0.0
    ξ₁::Float64 = 30.0
    scale::Float64 = 1.0 # Used for scaling to get better numerical stability
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

function abc(κ, ϵ, ω, λ)
    (; d, σ) = λ
    a = (1 / σ + im * ω / κ) / 2
    b = oftype(a, d) / 2
    c = -im * κ / (1 - im * ϵ) / 2
    return a, b, c
end

function P(ξ, κ, ϵ, ω, λ)
    a, b, c = abc(κ, ϵ, ω, λ)

    return U(a, b, c * ξ^2)
end

function P_dξ(ξ, κ, ϵ, ω, λ)
    a, b, c = abc(κ, ϵ, ω, λ)
    z_dξ = 2c * ξ

    return U_dz(a, b, c * ξ^2) * z_dξ
end

function E(ξ, κ, ϵ, ω, λ)
    a, b, c = abc(κ, ϵ, ω, λ)

    z = c * ξ^2

    return exp(z) * U(b - a, b, -z)
end

function E_dξ(ξ, κ, ϵ, ω, λ)
    a, b, c = abc(κ, ϵ, ω, λ)

    z = c * ξ^2
    z_dξ = 2c * ξ

    return exp(z) * (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ
end

function W(ξ, κ, ϵ, ω, λ)
    a, b, c = abc(κ, ϵ, ω, λ)

    z = c * ξ^2

    return -im * κ / (1 - im * ϵ) *
           exp(sign(imag(c)) * im * (b - a) * π) *
           ξ *
           z^-b *
           exp(z)
end

function B_W(κ, ϵ, ω, λ)
    (; δ) = λ

    a, b, c = abc(κ, ϵ, ω, λ)

    return -(1 + im * δ) / (im * κ) * exp(-sign(imag(c)) * im * (b - a) * π) * c^b
end

function J_P(ξ, κ, ϵ, ω, λ)
    (; δ) = λ

    return (1 + im * δ) / (1 - im * ϵ) * P(ξ, κ, ϵ, ω, λ) / W(ξ, κ, ϵ, ω, λ)
end

function J_E(ξ, κ, ϵ, ω, λ)
    (; δ) = λ

    return (1 + im * δ) / (1 - im * ϵ) * E(ξ, κ, ϵ, ω, λ) / W(ξ, κ, ϵ, ω, λ)
end

# Optimized for the case d == 1 and δ = 0
function system_d1_δ0(u, (κ, ϵ, ω, λ), ξ)
    (; d, σ) = λ
    a, b, α, β = u

    @fastmath begin
        a2b2σ = (a^2 + b^2)^σ

        F1 = κ * ξ * β + κ / σ * b + ω * a - a2b2σ * a
        F2 = -κ * ξ * α - κ / σ * a + ω * b - a2b2σ * b

        return SVector(α, β, (F1 - ϵ * F2) / (1 + ϵ^2), (ϵ * F1 + F2) / (1 + ϵ^2))
    end
end

# Optimized for the case d == 3, σ == 1 and δ = 0
function system_d3_σ1_δ0(u, (κ, ϵ, ω, λ), ξ)
    a, b, α, β = u

    @fastmath begin
        a2b2σ = a^2 + b^2

        F1 = κ * ξ * β + κ * b + ω * a - a2b2σ * a
        F2 = -κ * ξ * α - κ * a + ω * b - a2b2σ * b

        if !iszero(ξ)
            F1 -= 2(α + ϵ * β) / ξ
            F2 -= 2(β - ϵ * α) / ξ
        end

        return SVector(α, β, (F1 - ϵ * F2) / (1 + ϵ^2), (ϵ * F1 + F2) / (1 + ϵ^2))
    end
end

function system(u, (κ, ϵ, ω, λ), ξ)
    (; d, σ, δ) = λ
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

function G(μ, κ, ϵ, ω, λ::Params)
    (; d, σ, ξ₁, scale) = λ

    μ, κ, ω = scale^(1 / σ) * μ, scale^2 * κ, scale^2 * ω

    if d == 1 && λ.δ == 0
        prob = ODEProblem{false}(
            system_d1_δ0,
            SVector(μ, 0, 0, 0),
            (zero(ξ₁), ξ₁),
            (κ, ϵ, ω, λ),
        )
        sol = solve(
            prob,
            AutoVern7(Rodas5P()),
            abstol = 1e-9,
            reltol = 1e-9,
            maxiters = 4000,
            save_everystep = false,
            verbose = false,
        )
    elseif d == 3 && σ == 1 && λ.δ == 0
        prob = ODEProblem{false}(
            system_d3_σ1_δ0,
            SVector(μ, 0, 0, 0),
            (zero(ξ₁), ξ₁),
            (κ, ϵ, ω, λ),
        )
        sol = solve(
            prob,
            AutoVern7(Rodas5P()),
            abstol = 1e-9,
            reltol = 1e-9,
            maxiters = 8000,
            save_everystep = false,
            verbose = false,
        )
    else
        prob = ODEProblem{false}(system, SVector(μ, 0, 0, 0), (zero(ξ₁), ξ₁), (κ, ϵ, ω, λ))
        sol = solve(
            prob,
            AutoVern7(Rodas5P()),
            abstol = 1e-9,
            reltol = 1e-9,
            maxiters = 8000,
            save_everystep = false,
            verbose = false,
        )
    end

    a, b, α, β = sol[end]::SVector{4,promote_type(typeof(κ), typeof(ϵ))}

    order = 2 # Order of approximation to use

    if sol.retcode != ReturnCode.Success
        # IMPROVE: We don't really want to compute anything in this
        # case. But it is unclear exactly what we should return then.
        # For now we at least only use the first order approximation.
        order = 1
    end

    if κ < 0
        # IMPROVE: We don't really want to compute anything in this
        # case. But it is unclear exactly what we should return then.
        # For now we at least only use the first order approximation.
        order = 1
    end

    Q_0, dQ_0 = complex(a, b), complex(α, β)

    if order == 1 # First order approximation
        γ = Q_0 / P(ξ₁, κ, ϵ, ω, λ)

        dQ_inf = γ * P_dξ(ξ₁, κ, ϵ, ω, λ)
    elseif order == 2 # Second order approximation
        p = P(ξ₁, κ, ϵ, ω, λ)
        p_dξ = P_dξ(ξ₁, κ, ϵ, ω, λ)
        e = E(ξ₁, κ, ϵ, ω, λ)
        e_dξ = E_dξ(ξ₁, κ, ϵ, ω, λ)

        _, _, c = abc(κ, ϵ, ω, λ)

        I_P_witout_γ =
            B_W(κ, ϵ, ω, λ) * exp(-c * ξ₁^2) * p * ξ₁^(d - 2) * abs(p)^2σ * p / 2c

        γ = nlsolve([Q_0 / p], ftol = 1e-14) do γ
            Q_0 - (γ[1] * p + e * I_P_witout_γ * abs(γ[1])^2σ * γ[1])
        end.zero[1]

        # Compute derivative for solution at infinity with given γ
        I_P = I_P_witout_γ * abs(γ)^2σ * γ

        Q_inf = γ * p + e * I_P

        I_E_dξ = J_E(ξ₁, κ, ϵ, ω, λ) * abs(Q_inf)^2σ * Q_inf
        I_P_dξ = -J_P(ξ₁, κ, ϵ, ω, λ) * abs(Q_inf)^2σ * Q_inf

        dQ_inf = γ * p_dξ + p * I_E_dξ + e_dξ * I_P + e * I_P_dξ
    else
        error("unsupported order $order")
    end

    res = dQ_0 - dQ_inf

    return [real(res), imag(res)]
end

G(x, (ϵ, ω, λ)::@NamedTuple{ϵ::S, ω::T, λ::Params}) where {S,T} = G(x[1], x[2], ϵ, ω, λ)
G(x, (mκ, ω, λ)::@NamedTuple{mκ::S, ω::T, λ::Params}) where {S,T} =
    -G(x[1], -mκ, x[2], ω, λ)
G(x, (κ, ω, λ)::@NamedTuple{κ::S, ω::T, λ::Params}) where {S,T} = G(x[1], κ, x[2], ω, λ)
G(x, (μ, ϵ, λ)::@NamedTuple{μ::S, ϵ::T, λ::Params}) where {S,T} = G(μ, x[1], ϵ, x[2], λ)
G(x, (μ, mκ, λ)::@NamedTuple{μ::S, mκ::T, λ::Params}) where {S,T} = G(μ, -mκ, x[1], x[2], λ)

# Like G but uses the first two terms in the asymptotic expansion
function G_asym(μ, κ, ϵ, ω, λ::Params)
    (; d, σ, ξ₁, scale) = λ

    μ, κ, ω = scale^(1 / σ) * μ, scale^2 * κ, scale^2 * ω

    if d == 1 && λ.δ == 0
        prob = ODEProblem{false}(
            system_d1_δ0,
            SVector(μ, 0, 0, 0),
            (zero(ξ₁), ξ₁),
            (κ, ϵ, ω, λ),
        )
        sol = solve(
            prob,
            AutoVern7(Rodas5P()),
            abstol = 1e-9,
            reltol = 1e-9,
            maxiters = 4000,
            save_everystep = false,
            verbose = false,
        )
    elseif d == 3 && σ == 1 && λ.δ == 0
        prob = ODEProblem{false}(
            system_d3_σ1_δ0,
            SVector(μ, 0, 0, 0),
            (zero(ξ₁), ξ₁),
            (κ, ϵ, ω, λ),
        )
        sol = solve(
            prob,
            AutoVern7(Rodas5P()),
            abstol = 1e-9,
            reltol = 1e-9,
            maxiters = 8000,
            save_everystep = false,
            verbose = false,
        )
    else
        prob = ODEProblem{false}(system, SVector(μ, 0, 0, 0), (zero(ξ₁), ξ₁), (κ, ϵ, ω, λ))
        sol = solve(
            prob,
            AutoVern7(Rodas5P()),
            abstol = 1e-9,
            reltol = 1e-9,
            maxiters = 8000,
            save_everystep = false,
            verbose = false,
        )
    end

    a, b, α, β = sol[end]::SVector{4,promote_type(typeof(κ), typeof(ϵ))}

    order = 2 # Order of approximation to use

    if sol.retcode != ReturnCode.Success
        # IMPROVE: We don't really want to compute anything in this
        # case. But it is unclear exactly what we should return then.
        # For now we at least only use the first order approximation.
        order = 1
    end

    Q_0, dQ_0 = complex(a, b), complex(α, β)

    if order == 1 # First order approximation
        # We want this:
        # c_0 = Q_0 / ξ₁^-2a
        # dQ_inf = -2a * c_0 * ξ₁^(-2a - 1)
        # Which simplifies to this:
        dQ_inf = -2a * Q_0 / ξ₁
    elseif order == 2 # Second order approximation
        a, b, c = abc(κ, ϵ, ω, λ)

        c_0 = nlsolve([Q_0 / ξ₁^-2a], ftol = 1e-14) do c_0
            Q_0 - (
                c_0[1] * ξ₁^-2a +
                c_0[1] * ξ₁^(-2a - 2) * (2a * (2a + 1) * (1 - im * ϵ) + abs(c_0[1])^2σ) / (2im * κ)
            )
        end.zero[1]

        # Compute derivative for solution at infinity with given c_0
        dQ_inf =
            -2a * c_0 * ξ₁^(-2a - 1) +
            (-2a - 2) * c_0 * ξ₁^(-2a - 3) * (2a * (2a + 1) * (1 - im * ϵ) + abs(c_0)^2σ) /
            (2im * κ)
    else
        error("unsupported order $order")
    end

    res = dQ_0 - dQ_inf

    return [real(res), imag(res)]
end

G_asym(x, (ϵ, ω, λ)::@NamedTuple{ϵ::S, ω::T, λ::Params}) where {S,T} =
    G_asym(x[1], x[2], ϵ, ω, λ)
G_asym(x, (mκ, ω, λ)::@NamedTuple{mκ::S, ω::T, λ::Params}) where {S,T} =
    G_asym(x[1], -mκ, x[2], ω, λ)
G_asym(x, (μ, ϵ, λ)::@NamedTuple{μ::S, ϵ::T, λ::Params}) where {S,T} =
    G_asym(μ, x[1], ϵ, x[2], λ)
G_asym(x, (μ, mκ, λ)::@NamedTuple{μ::S, mκ::T, λ::Params}) where {S,T} =
    G_asym(μ, -mκ, x[1], x[2], λ)

function sverak_initial(j, d; fix_omega = true, autoscale = fix_omega)
    if d == 1
        μs = [1.23204, 0.78308, 1.12389, 0.88393, 1.07969, 0.92761, 1.05440, 0.94914]
        κs = [0.85310, 0.49322, 0.34680, 0.26678, 0.21621, 0.18192, 0.15667, 0.13749]
        λ = Params()
    elseif d == 3
        μs = [
            1.885903265965844,
            0.839882186873777,
            1.1174528182552068,
            0.9713988670269424,
            1.0210957154390576,
        ]
        κs = [
            0.9174206192134661,
            0.3212555262690968,
            0.22690292062401618,
            0.1698342794838976,
            0.13810720151446634,
        ]
        ξ₁ = j == 5 ? 50.0 : 30.0
        scale = autoscale ? [0.5, 0.8, 1.0, 1.0, 0.85][j] : 1
        λ = Params(d = 3, σ = 1; ξ₁, scale)
    end

    if fix_omega
        return μs[j], κs[j], 0.0, 1.0, λ
    else
        scaling = μs[j]^-λ.σ
        return 1.0, scaling^2 * κs[j], 0.0, scaling^2, λ
    end
end

# Bifurcation with ϵ as the bifurcation parameter
function branch_epsilon(
    μ,
    κ,
    ϵ,
    ω,
    λ::Params;
    asym = false,
    fix_omega = true,
    max_steps = nothing,
    ϵ_stop = nothing,
)
    # IMPROVE: Look at using ShootingProblem
    #prob = BifurcationProblem(
    #    ifelse(asym, G_asym, G),
    #    [μ, κ],
    #    (; ϵ, ω, λ),
    #    (@lens _.ϵ),
    #    record_from_solution = (x, _) -> (κ = x[2], μ = x[1]),
    #)

    if fix_omega
        prob = BifurcationProblem(
            ifelse(asym, G_asym, G),
            [μ, κ],
            (; ϵ, ω, λ),
            (@lens _.ϵ),
            record_from_solution = (x, _) -> (κ = x[2], μ = x[1]),
        )
    else
        prob = BifurcationProblem(
            ifelse(asym, G_asym, G),
            [κ, ω],
            (; μ, ϵ, λ),
            (@lens _.ϵ),
            record_from_solution = (x, _) -> (κ = x[1], ω = x[2]),
        )
    end

    if λ.d == 1
        opts = ContinuationPar(
            p_min = 0.0,
            p_max = 0.07,
            dsmin = 5e-5,
            ds = 1e-4,
            dsmax = 5e-4,
            max_steps = something(max_steps, 1500),
            detect_bifurcation = 0,
            newton_options = NewtonPar(tol = 1e-6, max_iterations = 10),
        )

        finalise_solution =
            (z, tau, step, contResult; kwargs...) ->
                !(tau.p < 0 && z.p < something(ϵ_stop, 0.02))
    elseif λ.d == 3
        # IMPROVE: This is not able to capture the full branches. It
        # needs more tuning or other changes.
        opts = ContinuationPar(
            p_min = 0.0,
            p_max = 0.3,
            dsmin = 5e-6,
            ds = 1e-4,
            dsmax = 5e-4,
            max_steps = something(max_steps, 15000),
            detect_bifurcation = 0,
            newton_options = NewtonPar(tol = 1e-6, max_iterations = 20),
        )

        finalise_solution =
            (z, tau, step, contResult; kwargs...) ->
                !(tau.p < 0 && z.p < something(ϵ_stop, 0.075))
    end

    # Parameters should always be positive
    callback_newton = (values; kwargs...) -> values.x[1] > 0 && values.x[2] > 0

    br = continuation(prob, PALC(), opts; finalise_solution, callback_newton)

    return br
end

# Bifurcation with κ as the bifurcation parameter
function branch_kappa(
    μ,
    κ,
    ϵ,
    ω,
    λ::Params;
    asym = false,
    fix_omega = true,
    max_steps = nothing,
    ϵ_stop = nothing,
)
    # IMPROVE: I haven't figure out how to chose the direction of the
    # initial bifurcation direction. To get it in the right direction
    # we work with -κ instead of κ.
    if fix_omega
        # For reasons I have not fully understood we need a minus sign
        # in front of G the non-asymptotic case to get the bifurcation
        # in the right direction.
        prob = BifurcationProblem(
            ifelse(asym, G_asym, (-) ∘ G),
            [μ, ϵ],
            (; mκ = -κ, ω, λ),
            (@lens _.mκ),
            record_from_solution = (x, _) -> (ϵ = x[2], μ = x[1]),
        )
    else
        prob = BifurcationProblem(
            ifelse(asym, G_asym, G),
            [ϵ, ω],
            (; μ, mκ = -κ, λ),
            (@lens _.mκ),
            record_from_solution = (x, _) -> (ϵ = x[1], ω = x[2]),
        )
    end

    if λ.d == 1
        opts = ContinuationPar(
            dsmin = 5e-5,
            ds = 1e-4,
            dsmax = 5e-4,
            max_steps = something(max_steps, 1500),
            detect_bifurcation = 0,
            newton_options = NewtonPar(tol = 1e-6, max_iterations = 10),
        )

        finalise_solution = if fix_omega
            (z, tau, step, contResult; kwargs...) ->
                !(tau.u[2] < 0 && z.u[2] < something(ϵ_stop, 0.02))
        else
            (z, tau, step, contResult; kwargs...) ->
                !(tau.u[1] < 0 && z.u[1] < something(ϵ_stop, 0.02))
        end
    elseif λ.d == 3
        opts = ContinuationPar(
            dsmin = 5e-6,
            ds = 1e-4,
            dsmax = 5e-4,
            max_steps = something(max_steps, 15000),
            detect_bifurcation = 0,
            newton_options = NewtonPar(tol = 1e-6, max_iterations = 10),
        )

        finalise_solution = if fix_omega
            (z, tau, step, contResult; kwargs...) ->
                !(tau.u[2] < 0 && z.u[2] < something(ϵ_stop, 0.075))
        else
            (z, tau, step, contResult; kwargs...) ->
                !(tau.u[1] < 0 && z.u[1] < something(ϵ_stop, 0.075))
        end
    end

    br = continuation(prob, PALC(), opts; finalise_solution)

    return br
end

end # module CGLBranch
