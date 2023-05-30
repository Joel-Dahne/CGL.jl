"""
    approximate_parameters_F(p::AbstractGLParams{T}, ξ₁::T) where {T}

Return a function `F(κ, μ)` for computing
```
u_dξ(ξ₁, (p, κ)) - γ * P_dξ(ξ₁, (p, κ))
```
where `γ = u(ξ₁, (p, κ)) / P(ξ₁, (p, κ))`.

Here `u(ξ, (p, κ))` is the solution satisfying the initial condition
`u(0, (p, κ)) = μ` and `u_dξ` is the derivative with respect to `ξ`.
"""
function approximate_parameters_F(p::AbstractGLParams{T}, ξ₁::T) where {T}
    return (κ, μ) -> begin
        prob = ODEProblem(gl_equation_real, SVector(μ, 0.0, 0.0, 0.0), (zero(ξ₁), ξ₁), (p, κ))

        u_solution = solve(prob, abstol = 1e-9, reltol = 1e-9)

        u = u_solution[end][1] + im * u_solution[end][2]
        u_dξ = u_solution[end][3] + im * u_solution[end][4]

        γ = u / P(ξ₁, (p, κ))

        return u_dξ - γ * P_dξ(ξ₁, (p, κ))
    end
end

"""
    approximate_parameters(κ₀::T, μ₀::T, p::AbstractGLParams{T}, ξ₁::T) where {T}

Compute `κ` and `μ` such that
```
u_dξ(ξ₁, (p, κ)) = γ * P_dξ(ξ₁, (p, κ))
```
where `γ = u(ξ₁, (p, κ)) / P(ξ₁, (p, κ))`.

Here `u(ξ, (p, κ))` is the solution satisfying the initial condition
`u(0, (p, κ)) = μ` and `u_dξ` is the derivative with respect to `ξ`.
"""
function approximate_parameters(κ₀::T, μ₀::T, p::AbstractGLParams{T}, ξ₁::T) where {T}
    F = splat(approximate_parameters_F(p, ξ₁))
    sol = nlsolve([κ₀, μ₀], xtol = 1e-12, ftol = 1e-12) do inp
        res = F(inp)
        return [real(res), imag(res)]
    end

    return sol.zero
end
