"""
    Q_zero(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Q_zero(μ::Arb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}; tol::Float64 = 1e-11)
    sol = Q_zero_capd(μ, κ, ϵ, ξ₁, λ; tol)
    return SVector(Acb(sol[1], sol[2]), Acb(sol[3], sol[4]))
end

function Q_zero(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    sol = Q_zero_float(μ, κ, ϵ, ξ₁, λ; tol)
    return SVector(complex(sol[1], sol[2]), complex(sol[3], sol[4]))
end

"""
    Q_zero_jacobian_kappa(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_zero`](@ref) w.r.t. the
parameters `μ` and `κ`.
"""
function Q_zero_jacobian_kappa(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    J = Q_zero_jacobian_kappa_capd(μ, κ, ϵ, ξ₁, λ; tol)

    return SMatrix{2,2}(
        Acb(J[1, 1], J[2, 1]),
        Acb(J[3, 1], J[4, 1]),
        Acb(J[1, 2], J[2, 2]),
        Acb(J[3, 2], J[4, 2]),
    )
end

function Q_zero_jacobian_kappa(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    J = Q_zero_jacobian_kappa_float(μ, κ, ϵ, ξ₁, λ; tol)

    return SMatrix{2,2}(
        complex(J[1, 1], J[2, 1]),
        complex(J[3, 1], J[4, 1]),
        complex(J[1, 2], J[2, 2]),
        complex(J[3, 2], J[4, 2]),
    )
end

"""
    Q_zero_jacobian_epsilon(μ, κ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_zero`](@ref) w.r.t. the
parameters `μ` and `ϵ`.
"""
function Q_zero_jacobian_epsilon(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    J = Q_zero_jacobian_epsilon_capd(μ, κ, ϵ, ξ₁, λ; tol)

    return SMatrix{2,2}(
        Acb(J[1, 1], J[2, 1]),
        Acb(J[3, 1], J[4, 1]),
        Acb(J[1, 2], J[2, 2]),
        Acb(J[3, 2], J[4, 2]),
    )
end

function Q_zero_jacobian_epsilon(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    J = Q_zero_jacobian_epsilon_float(μ, κ, ϵ, ξ₁, λ; tol)

    return SMatrix{2,2}(
        complex(J[1, 1], J[2, 1]),
        complex(J[3, 1], J[4, 1]),
        complex(J[1, 2], J[2, 2]),
        complex(J[3, 2], J[4, 2]),
    )
end
