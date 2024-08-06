"""
    Q_zero_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with four real values, the first two are the real and imaginary
values at `ξ₁` and the second two are their derivatives.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `Float64`. However, for `Arb` input with wide
intervals for `μ`, `κ` and/or `ϵ` it computes it at the corners of the
box they form. This means you still get something that resembles an
enclosure.
"""
function Q_zero_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    prob =
        ODEProblem{false}(cgl_equation_real, SVector(μ, 0, 0, 0), (zero(ξ₁), ξ₁), (κ, ϵ, λ))

    # Used to exit early in extreme cases. Sometimes the
    # refine_approximation methods en up in a bad region and the
    # solutions are highly oscillating, which makes the solver
    # extremely slow. For that reason we exit early in case the norm
    # of the solution grows to large.
    unstable_check = (dt, u, p, t) -> any(isnan, u) || norm(u) > 100

    sol = solve(
        prob,
        AutoVern7(Rodas5P()),
        abstol = tol,
        reltol = tol,
        save_everystep = false,
        verbose = false;
        unstable_check,
    )

    return sol.u[end]
end

function Q_zero_float(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    ξ₁ = Float64(ξ₁)
    λ = CGLParams{Float64}(λ)

    μs = iswide(μ) ? collect(Float64.(getinterval(μ))) : [Float64(μ)]
    κs = iswide(κ) ? collect(Float64.(getinterval(κ))) : [Float64(κ)]
    ϵs = iswide(ϵ) ? collect(Float64.(getinterval(ϵ))) : [Float64(ϵ)]

    us = map(Iterators.product(μs, κs, ϵs)) do (μ, κ, ϵ)
        Q_zero_float(μ, κ, ϵ, ξ₁, λ; tol)
    end

    return SVector(
        Arb(extrema(getindex.(us, 1))),
        Arb(extrema(getindex.(us, 2))),
        Arb(extrema(getindex.(us, 3))),
        Arb(extrema(getindex.(us, 4))),
    )
end

"""
    Q_zero_jacobian_kappa_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_zero_float`](@ref) w.r.t.
the parameters `μ` and `κ`.
"""
function Q_zero_jacobian_kappa_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    return ForwardDiff.jacobian(SVector(μ, κ)) do (μ, κ)
        Q_zero_float(μ, κ, ϵ, ξ₁, λ; tol)
    end
end

function Q_zero_jacobian_kappa_float(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    ξ₁ = Float64(ξ₁)
    λ = CGLParams{Float64}(λ)

    μs = iswide(μ) ? collect(Float64.(getinterval(μ))) : [Float64(μ)]
    κs = iswide(κ) ? collect(Float64.(getinterval(κ))) : [Float64(κ)]
    ϵs = iswide(ϵ) ? collect(Float64.(getinterval(ϵ))) : [Float64(ϵ)]

    Js = map(Iterators.product(μs, κs, ϵs)) do (μ, κ, ϵ)
        Q_zero_jacobian_kappa_float(μ, κ, ϵ, ξ₁, λ; tol)
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
    Q_zero_jacobian_epsilon_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_zero_float`](@ref) w.r.t.
the parameters `μ` and `ϵ`.
"""
function Q_zero_jacobian_epsilon_float(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    return ForwardDiff.jacobian(SVector(μ, ϵ)) do (μ, ϵ)
        Q_zero_float(μ, κ, ϵ, ξ₁, λ; tol)
    end
end

function Q_zero_jacobian_epsilon_float(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    ξ₁ = Float64(ξ₁)
    λ = CGLParams{Float64}(λ)

    μs = iswide(μ) ? collect(Float64.(getinterval(μ))) : [Float64(μ)]
    κs = iswide(κ) ? collect(Float64.(getinterval(κ))) : [Float64(κ)]
    ϵs = iswide(ϵ) ? collect(Float64.(getinterval(ϵ))) : [Float64(ϵ)]

    Js = map(Iterators.product(μs, κs, ϵs)) do (μ, κ, ϵ)
        Q_zero_jacobian_epsilon_float(μ, κ, ϵ, ξ₁, λ; tol)
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
    Q_zero_float_curve(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Similar to [`Q_zero_float`](@ref) but returns the whole
solution object given by the ODE solver, instead of just the value at
the final point.
"""
function Q_zero_float_curve(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    prob =
        ODEProblem{false}(cgl_equation_real, SVector(μ, 0, 0, 0), (zero(ξ₁), ξ₁), (κ, ϵ, λ))

    sol = solve(prob, AutoVern7(Rodas5P()), abstol = tol, reltol = tol, verbose = false)

    return sol
end
