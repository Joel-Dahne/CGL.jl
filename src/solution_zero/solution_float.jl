"""
    solution_zero_float(μ::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes `u(ξ₁)`.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `Float64`. However, for `T = Arb` with wide
intervals for `μ` and/or `κ` it computes it at the corners of the box
they form. This means you still get something that resembles an
enclosure.
"""
function solution_zero_float(μ, κ, ξ₁, λ::CGLParams)
    prob = ODEProblem{false}(
        cgl_equation_real_alt,
        SVector(μ, 0, 0, 0),
        (zero(ξ₁), ξ₁),
        (κ, λ),
    )

    sol = solve(
        prob,
        AutoVern7(Rodas5P()),
        abstol = 1e-10,
        reltol = 1e-10,
        save_everystep = false,
        verbose = false,
    )

    return sol.u[end]
end

function solution_zero_float(μ::Arb, κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
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

    us = map(Iterators.product(μs, κs)) do (μ, κ)
        solution_zero_float(μ, κ, ξ₁, λ)
    end

    return SVector(
        Arb(extrema(getindex.(us, 1))),
        Arb(extrema(getindex.(us, 2))),
        Arb(extrema(getindex.(us, 3))),
        Arb(extrema(getindex.(us, 4))),
    )
end

"""
    solution_zero_jacobian_float(μ::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes `u(ξ₁)` as well as the Jacobian w.r.t. `μ` and
`κ`.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `Float64`. However, for `T = Arb` with wide
intervals for `μ` and/or `κ` it computes it at the corners of the box
they form. This means you still get something that resembles an
enclosure.
"""
function solution_zero_jacobian_float(
    μ::Float64,
    κ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    u = solution_zero_float(μ, κ, ξ₁, λ)

    J = ForwardDiff.jacobian(SVector(μ, κ)) do (μ, κ)
        solution_zero_float(μ, κ, ξ₁, λ)
    end

    return u, J
end

function solution_zero_jacobian_float(μ::Arb, κ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
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

    res = map(Iterators.product(μs, κs)) do (μ, κ)
        solution_zero_jacobian_float(μ, κ, ξ₁, λ)
    end

    us = getindex.(res, 1)
    jacobians = getindex.(res, 2)

    u = SVector(
        Arb(extrema(getindex.(us, 1))),
        Arb(extrema(getindex.(us, 2))),
        Arb(extrema(getindex.(us, 3))),
        Arb(extrema(getindex.(us, 4))),
    )

    J = SMatrix{4,2}(
        Arb(extrema(getindex.(jacobians, 1))),
        Arb(extrema(getindex.(jacobians, 2))),
        Arb(extrema(getindex.(jacobians, 3))),
        Arb(extrema(getindex.(jacobians, 4))),
        Arb(extrema(getindex.(jacobians, 5))),
        Arb(extrema(getindex.(jacobians, 6))),
        Arb(extrema(getindex.(jacobians, 7))),
        Arb(extrema(getindex.(jacobians, 8))),
    )

    return u, J
end
