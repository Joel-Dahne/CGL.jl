"""
    G(μ::T, γ_real::T, γ_imag::T, κ::T, ϵ::T, ξ₁::T, λ::CGLParams{T}) where {T}

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
function G(μ::T, γ_real::T, γ_imag::T, κ::T, ϵ::T, ξ₁::T, λ::CGLParams{T}) where {T}
    Q_0, Q_0_dξ = solution_zero(μ, κ, ϵ, ξ₁, λ)

    γ = T == Arb ? Acb(γ_real, γ_imag) : complex(γ_real, γ_imag)
    Q_inf, Q_inf_dξ = solution_infinity(γ, κ, ϵ, ξ₁, λ)

    G1 = Q_0 - Q_inf
    G2 = Q_0_dξ - Q_inf_dξ

    SVector(real(G1), imag(G1), real(G2), imag(G2))
end

"""
    G_jacobian_kappa(μ::T, γ_real::T, γ_imag::T, κ::T, ϵ::T, ξ₁::T, λ::CGLParams{T}) where {T}

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
function G_jacobian_kappa(
    μ::T,
    γ_real::T,
    γ_imag::T,
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    # TODO: Allowing more control of when to use the mincing version
    # and whether it uses threading or not.
    if λ.d == 3 && iszero(ϵ) && (iswide(μ) || iswide(κ))
        μs = mince(μ, ifelse(iswide(μ), 4, 1))
        κs = mince(κ, ifelse(iswide(κ), 96, 1))

        # IMPROVE: We currently fix a tighter tolerance here. This
        # should possibly be adjustable.
        Q_0_Js = tmap(
            ((μ, κ),) -> solution_zero_jacobian_kappa(μ, κ, ϵ, ξ₁, λ, tol = 1e-14),
            collect(Iterators.product(μs, κs)),
        )
        Q_0_J = SMatrix{2,2}(
            Arblib.union(getindex.(Q_0_Js, 1)...),
            Arblib.union(getindex.(Q_0_Js, 2)...),
            Arblib.union(getindex.(Q_0_Js, 3)...),
            Arblib.union(getindex.(Q_0_Js, 4)...),
        )

        γ = T == Arb ? Acb(γ_real, γ_imag) : complex(γ_real, γ_imag)

        Q_inf_Js = tmap(κ -> solution_infinity_jacobian_kappa(γ, κ, ϵ, ξ₁, λ), κs)
        Q_inf_J = SMatrix{2,2}(
            Arblib.union(getindex.(Q_inf_Js, 1)...),
            Arblib.union(getindex.(Q_inf_Js, 2)...),
            Arblib.union(getindex.(Q_inf_Js, 3)...),
            Arblib.union(getindex.(Q_inf_Js, 4)...),
        )
    else
        Q_0_J = solution_zero_jacobian_kappa(μ, κ, ϵ, ξ₁, λ)

        γ = T == Arb ? Acb(γ_real, γ_imag) : complex(γ_real, γ_imag)
        Q_inf_J = solution_infinity_jacobian_kappa(γ, κ, ϵ, ξ₁, λ)
    end

    return SMatrix{4,4,T}(
        real(Q_0_J[1, 1]),
        imag(Q_0_J[1, 1]),
        real(Q_0_J[2, 1]),
        imag(Q_0_J[2, 1]),
        # Derivatives w.r.t γ_real
        -real(Q_inf_J[1, 1]),
        -imag(Q_inf_J[1, 1]),
        -real(Q_inf_J[2, 1]),
        -imag(Q_inf_J[2, 1]),
        # Derivatives w.r.t γ_imag
        imag(Q_inf_J[1, 1]),
        -real(Q_inf_J[1, 1]),
        imag(Q_inf_J[2, 1]),
        -real(Q_inf_J[2, 1]),
        # Derivatives w.r.t κ
        real(Q_0_J[1, 2]) - real(Q_inf_J[1, 2]),
        imag(Q_0_J[1, 2]) - imag(Q_inf_J[1, 2]),
        real(Q_0_J[2, 2]) - real(Q_inf_J[2, 2]),
        imag(Q_0_J[2, 2]) - imag(Q_inf_J[2, 2]),
    )
end

"""
    G_jacobian_epsilon(μ::T, γ_real::T, γ_imag::T, κ::T, ϵ::T, ξ₁::T, λ::CGLParams{T}) where {T}

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
function G_jacobian_epsilon(
    μ::T,
    γ_real::T,
    γ_imag::T,
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    Q_0_J = solution_zero_jacobian_epsilon(μ, κ, ϵ, ξ₁, λ)

    γ = T == Arb ? Acb(γ_real, γ_imag) : complex(γ_real, γ_imag)
    Q_inf_J = solution_infinity_jacobian_epsilon(γ, κ, ϵ, ξ₁, λ)

    return SMatrix{4,4,T}(
        real(Q_0_J[1, 1]),
        imag(Q_0_J[1, 1]),
        real(Q_0_J[2, 1]),
        imag(Q_0_J[2, 1]),
        # Derivatives w.r.t γ_real
        -real(Q_inf_J[1, 1]),
        -imag(Q_inf_J[1, 1]),
        -real(Q_inf_J[2, 1]),
        -imag(Q_inf_J[2, 1]),
        # Derivatives w.r.t γ_imag
        imag(Q_inf_J[1, 1]),
        -real(Q_inf_J[1, 1]),
        imag(Q_inf_J[2, 1]),
        -real(Q_inf_J[2, 1]),
        # Derivatives w.r.t κ
        real(Q_0_J[1, 2]) - real(Q_inf_J[1, 2]),
        imag(Q_0_J[1, 2]) - imag(Q_inf_J[1, 2]),
        real(Q_0_J[2, 2]) - real(Q_inf_J[2, 2]),
        imag(Q_0_J[2, 2]) - imag(Q_inf_J[2, 2]),
    )
end
