"""
    G_real(μ::T, γ_real::T, γ_imag::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}

Let `Q_0` be the solution to [`ivp_zero_complex`](@ref) and `Q_inf` a
solution to [`fpp_infinity_complex`](@ref), with `γ = γ_real + im *
γ_imag`. Let
```
G(ξ) = Q_0(ξ) - Q_inf(ξ)
```
This function computes `G(ξ₁)` and `d(G)(ξ₁)` in real variables, that
is
```
[real(G(ξ₁)), imag(G(ξ₁)), real(d(G)(ξ₁)), imag(d(G)(ξ₁))]
```
We here use `d(G)` to denote the derivative of `G` with respect to
`ξ`.
"""
function G_real(μ::T, γ_real::T, γ_imag::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}
    Q_0, Q_0_dξ = solution_zero(μ, κ, ξ₁, λ)

    γ = if T == Arb
        Acb(γ_real, γ_imag)
    else
        γ_real + im * γ_imag
    end
    Q_inf, Q_inf_dξ = solution_infinity(γ, κ, ξ₁, λ)

    G1 = Q_0 - Q_inf
    G2 = Q_0_dξ - Q_inf_dξ

    SVector(real(G1), imag(G1), real(G2), imag(G2))
end

"""
    G_jacobian_real(μ::T, γ_real::T, γ_imag::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}

Let `Q_0` be the solution to [`ivp_zero_complex`](@ref) and `Q_inf` a
solution to [`fpp_infinity_complex`](@ref), with `γ = γ_real + im *
γ_imag`. Let
```
G(ξ) = Q_0(ξ) - Q_inf(ξ)
```
This function computes the Jacobian of `G(ξ₁)` and `d(G)(ξ₁)` with
respect to `μ, γ_real, γ_imag, κ` in real variables. That is
```
[
    real(d(G(ξ₁), μ))    real(d(G(ξ₁), γ_)(ξ₁))    real(d(G(ξ₁), γ_)(ξ₁))    real(d(G(ξ₁), κ))
    imag(d(G(ξ₁), μ))    imag(d(G(ξ₁), γ_)(ξ₁))    imag(d(G(ξ₁), γ_)(ξ₁))    imag(d(G(ξ₁), κ))
    real(d(d(G)(ξ₁), μ)) real(d(d(G)(ξ₁), γ_)(ξ₁)) real(d(d(G)(ξ₁), γ_)(ξ₁)) real(d(d(G)(ξ₁), κ))
    imag(d(d(G)(ξ₁), μ)) imag(d(d(G)(ξ₁), γ_)(ξ₁)) imag(d(d(G)(ξ₁), γ_)(ξ₁)) imag(d(d(G)(ξ₁), κ))
]
```
We here use `d(G)` to denote the derivative of `G` with respect to `ξ`
and `d(G(ξ₁), μ)` to denote the derivative of `G(ξ₁)` with respect to
`μ`.
"""
function G_jacobian_real(μ::T, γ_real::T, γ_imag::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}
    Q_0_J = solution_zero_jacobian(μ, κ, ξ₁, λ)

    γ = if T == Arb
        Acb(γ_real, γ_imag)
    else
        γ_real + im * γ_imag
    end
    Q_inf_J = solution_infinity_jacobian(γ, κ, ξ₁, λ)

    G_J = zeros(T, 4, 4)

    # Derivatives w.r.t μ
    # Q_inf doesn't depend on μ so that derivative is zero
    G_J[1, 1] = real(Q_0_J[1, 1])
    G_J[2, 1] = imag(Q_0_J[1, 1])
    G_J[3, 1] = real(Q_0_J[2, 1])
    G_J[4, 1] = imag(Q_0_J[2, 1])

    # Derivatives w.r.t γ_real
    G_J[1, 2] = -real(Q_inf_J[1, 1])
    G_J[2, 2] = -imag(Q_inf_J[1, 1])
    G_J[3, 2] = -real(Q_inf_J[2, 1])
    G_J[4, 2] = -imag(Q_inf_J[2, 1])

    # Derivatives w.r.t γ_imag
    G_J[1, 3] = imag(Q_inf_J[1, 1])
    G_J[2, 3] = -real(Q_inf_J[1, 1])
    G_J[3, 3] = imag(Q_inf_J[2, 1])
    G_J[4, 3] = -real(Q_inf_J[2, 1])

    # Derivatives w.r.t κ
    G_J[1, 4] = real(Q_0_J[1, 2]) - real(Q_inf_J[1, 2])
    G_J[2, 4] = imag(Q_0_J[1, 2]) - imag(Q_inf_J[1, 2])
    G_J[3, 4] = real(Q_0_J[2, 2]) - real(Q_inf_J[2, 2])
    G_J[4, 4] = imag(Q_0_J[2, 2]) - imag(Q_inf_J[2, 2])

    return G_J
end

"""
    G_jacobian_epsilon_real(μ::T, γ_real::T, γ_imag::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}

Let `Q_0` be the solution to [`ivp_zero_complex`](@ref) and `Q_inf` a
solution to [`fpp_infinity_complex`](@ref), with `γ = γ_real + im *
γ_imag`. Let
```
G(ξ) = Q_0(ξ) - Q_inf(ξ)
```
This function computes the Jacobian of `G(ξ₁)` and `d(G)(ξ₁)` with
respect to `μ, γ_real, γ_imag, ϵ` in real variables. That is
```
[
  real(d(G(ξ₁), μ))    real(d(G(ξ₁), γ_)(ξ₁))    real(d(G(ξ₁), γ_)(ξ₁))    real(d(G(ξ₁), ϵ))
  imag(d(G(ξ₁), μ))    imag(d(G(ξ₁), γ_)(ξ₁))    imag(d(G(ξ₁), γ_)(ξ₁))    imag(d(G(ξ₁), ϵ))
  real(d(d(G)(ξ₁), μ)) real(d(d(G)(ξ₁), γ_)(ξ₁)) real(d(d(G)(ξ₁), γ_)(ξ₁)) real(d(d(G)(ξ₁), ϵ))
  imag(d(d(G)(ξ₁), μ)) imag(d(d(G)(ξ₁), γ_)(ξ₁)) imag(d(d(G)(ξ₁), γ_)(ξ₁)) imag(d(d(G)(ξ₁), ϵ))
]
```
We here use `d(G)` to denote the derivative of `G` with respect to `ξ`
and `d(G(ξ₁), μ)` to denote the derivative of `G(ξ₁)` with respect to
`μ`.
"""
function G_jacobian_epsilon_real(
    μ::T,
    γ_real::T,
    γ_imag::T,
    κ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    Q_0_J = solution_zero_jacobian_epsilon(μ, κ, ξ₁, λ)

    γ = if T == Arb
        Acb(γ_real, γ_imag)
    else
        γ_real + im * γ_imag
    end
    Q_inf_J = solution_infinity_jacobian_epsilon(γ, κ, ξ₁, λ)

    G_J = zeros(T, 4, 4)

    # Derivatives w.r.t μ
    # Q_inf doesn't depend on μ so that derivative is zero
    G_J[1, 1] = real(Q_0_J[1, 1])
    G_J[2, 1] = imag(Q_0_J[1, 1])
    G_J[3, 1] = real(Q_0_J[2, 1])
    G_J[4, 1] = imag(Q_0_J[2, 1])

    # Derivatives w.r.t γ_real
    G_J[1, 2] = -real(Q_inf_J[1, 1])
    G_J[2, 2] = -imag(Q_inf_J[1, 1])
    G_J[3, 2] = -real(Q_inf_J[2, 1])
    G_J[4, 2] = -imag(Q_inf_J[2, 1])

    # Derivatives w.r.t γ_imag
    G_J[1, 3] = imag(Q_inf_J[1, 1])
    G_J[2, 3] = -real(Q_inf_J[1, 1])
    G_J[3, 3] = imag(Q_inf_J[2, 1])
    G_J[4, 3] = -real(Q_inf_J[2, 1])

    # Derivatives w.r.t ϵ
    G_J[1, 4] = real(Q_0_J[1, 2]) - real(Q_inf_J[1, 2])
    G_J[2, 4] = imag(Q_0_J[1, 2]) - imag(Q_inf_J[1, 2])
    G_J[3, 4] = real(Q_0_J[2, 2]) - real(Q_inf_J[2, 2])
    G_J[4, 4] = imag(Q_0_J[2, 2]) - imag(Q_inf_J[2, 2])

    return G_J
end
