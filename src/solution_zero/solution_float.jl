"""
    solution_zero_float_solution(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Similar to [`solution_zero_float`](@ref) but returns the whole
solution object given by the ODE solver, instead of just the value at
the final point.
"""
function solution_zero_float_solution(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    prob = ODEProblem{false}(
        cgl_equation_real_alt,
        SVector(μ, 0, 0, 0),
        (zero(ξ₁), ξ₁),
        (κ, ϵ, λ),
    )

    sol = solve(prob, AutoVern7(Rodas5P()), abstol = tol, reltol = tol, verbose = false)

    return sol
end

"""
    solution_zero_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes `u(ξ₁)`.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `Float64`. However, for `Arb` input with wide
intervals for `μ`, `κ` and/or `ϵ` it computes it at the corners of the
box they form. This means you still get something that resembles an
enclosure.
"""
function solution_zero_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    prob = ODEProblem{false}(
        cgl_equation_real_alt,
        SVector(μ, 0, 0, 0),
        (zero(ξ₁), ξ₁),
        (κ, ϵ, λ),
    )

    sol = solve(
        prob,
        AutoVern7(Rodas5P()),
        abstol = tol,
        reltol = tol,
        save_everystep = false,
        verbose = false,
    )

    return sol.u[end]
end

function solution_zero_float(μ::Arb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}; tol)
    ξ₁ = Float64(ξ₁)
    λ = CGLParams{Float64}(λ)

    if iswide(μ)
        μs = collect(Float64.(getinterval(μ)))
    else
        μs = [Float64(μ)]
    end
    if iswide(κ)
        κs = collect(Float64.(getinterval(κ)))
    else
        κs = [Float64(κ)]
    end
    if iswide(ϵ)
        ϵs = collect(Float64.(getinterval(ϵ)))
    else
        ϵs = [Float64(ϵ)]
    end

    us = map(Iterators.product(μs, κs, ϵs)) do (μ, κ, ϵ)
        solution_zero_float(μ, κ, ϵ, ξ₁, λ; tol)
    end

    return SVector(
        Arb(extrema(getindex.(us, 1))),
        Arb(extrema(getindex.(us, 2))),
        Arb(extrema(getindex.(us, 3))),
        Arb(extrema(getindex.(us, 4))),
    )
end

"""
    solution_zero_jacobian_kappa_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes the Jacobian at `ξ₁` w.r.t. `μ` and `κ`.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `Float64`. However, for `Arb` input with wide
intervals for `μ`, `κ` and/or `ϵ` it computes it at the corners of the
box they form. This means you still get something that resembles an
enclosure.
"""
function solution_zero_jacobian_kappa_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    return ForwardDiff.jacobian(SVector(μ, κ)) do (μ, κ)
        solution_zero_float(μ, κ, ϵ, ξ₁, λ; tol)
    end
end

function solution_zero_jacobian_kappa_float(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    ξ₁ = Float64(ξ₁)
    λ = CGLParams{Float64}(λ)

    if iswide(μ)
        μs = collect(Float64.(getinterval(μ)))
    else
        μs = [Float64(μ)]
    end
    if iswide(κ)
        κs = collect(Float64.(getinterval(κ)))
    else
        κs = [Float64(κ)]
    end
    if iswide(ϵ)
        ϵs = collect(Float64.(getinterval(ϵ)))
    else
        ϵs = [Float64(ϵ)]
    end

    Js = map(Iterators.product(μs, κs, ϵs)) do (μ, κ, ϵ)
        solution_zero_jacobian_kappa_float(μ, κ, ϵ, ξ₁, λ; tol)
    end

    return SMatrix{4,2}(
        Arb(extrema(getindex.(Js, 1))),
        Arb(extrema(getindex.(Js, 2))),
        Arb(extrema(getindex.(Js, 3))),
        Arb(extrema(getindex.(Js, 4))),
        Arb(extrema(getindex.(Js, 5))),
        Arb(extrema(getindex.(Js, 6))),
        Arb(extrema(getindex.(Js, 7))),
        Arb(extrema(getindex.(Js, 8))),
    )
end

"""
    solution_zero_jacobian_epsilon_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes the Jacobian at `ξ₁` w.r.t. `μ` and `ϵ`.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `Float64`. However, for `Arb` input with wide
intervals for `μ`, `κ` and/or `ϵ` it computes it at the corners of the
box they form. This means you still get something that resembles an
enclosure.
"""
function solution_zero_jacobian_epsilon_float(
    μ,
    κ,
    ϵ,
    ξ₁,
    λ::CGLParams;
    tol::Float64 = 1e-11,
)
    return ForwardDiff.jacobian(SVector(μ, ϵ)) do (μ, ϵ)
        solution_zero_float(μ, κ, ϵ, ξ₁, λ; tol)
    end
end

function solution_zero_jacobian_epsilon_float(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    ξ₁ = Float64(ξ₁)
    λ = CGLParams{Float64}(λ)

    if iswide(μ)
        μs = collect(Float64.(getinterval(μ)))
    else
        μs = [Float64(μ)]
    end
    if iswide(κ)
        κs = collect(Float64.(getinterval(κ)))
    else
        κs = [Float64(κ)]
    end
    if iswide(ϵ)
        ϵs = collect(Float64.(getinterval(ϵ)))
    else
        ϵs = [Float64(ϵ)]
    end

    Js = map(Iterators.product(μs, κs, ϵs)) do (μ, κ, ϵ)
        solution_zero_jacobian_epsilon_float(μ, κ, ϵ, ξ₁, λ; tol)
    end

    return SMatrix{4,2}(
        Arb(extrema(getindex.(Js, 1))),
        Arb(extrema(getindex.(Js, 2))),
        Arb(extrema(getindex.(Js, 3))),
        Arb(extrema(getindex.(Js, 4))),
        Arb(extrema(getindex.(Js, 5))),
        Arb(extrema(getindex.(Js, 6))),
        Arb(extrema(getindex.(Js, 7))),
        Arb(extrema(getindex.(Js, 8))),
    )
end
