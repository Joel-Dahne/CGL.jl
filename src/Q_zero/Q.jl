"""
    Q_zero(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Q_zero(μ::T, κ::T, ϵ::T, ξ₁::T, λ::CGLParams{T}; tol::Float64 = 1e-11) where {T}
    Q = if T == Arb
        Q_zero_capd(μ, κ, ϵ, ξ₁, λ; tol)
    else
        Q_zero_float(μ, κ, ϵ, ξ₁, λ; tol)
    end
    return SVector(_complex(Q[1], Q[2]), _complex(Q[3], Q[4]))
end

"""
    Q_zero_jacobian_kappa(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_zero`](@ref) w.r.t. the
parameters `μ` and `κ`.
"""
function Q_zero_jacobian_kappa(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    J = if T == Arb
        Q_zero_jacobian_kappa_capd(μ, κ, ϵ, ξ₁, λ; tol)
    else
        Q_zero_jacobian_kappa_float(μ, κ, ϵ, ξ₁, λ; tol)
    end

    return SMatrix{2,2}(
        _complex(J[1, 1], J[2, 1]),
        _complex(J[3, 1], J[4, 1]),
        _complex(J[1, 2], J[2, 2]),
        _complex(J[3, 2], J[4, 2]),
    )
end

"""
    Q_zero_jacobian_epsilon(μ, κ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_zero`](@ref) w.r.t. the
parameters `μ` and `ϵ`.
"""
function Q_zero_jacobian_epsilon(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    J = if T == Arb
        Q_zero_jacobian_epsilon_capd(μ, κ, ϵ, ξ₁, λ; tol)
    else
        Q_zero_jacobian_epsilon_float(μ, κ, ϵ, ξ₁, λ; tol)
    end

    return SMatrix{2,2}(
        _complex(J[1, 1], J[2, 1]),
        _complex(J[3, 1], J[4, 1]),
        _complex(J[1, 2], J[2, 2]),
        _complex(J[3, 2], J[4, 2]),
    )
end
