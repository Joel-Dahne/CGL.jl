"""
    G(μ, γ_real, γ_imag, κ, ϵ, ξ₁, λ::CGLParams)

Let
```
G(ξ) = Q_0(ξ) - Q_inf(ξ)
```
where `Q_0` is given by [`Q_zero`](@ref) and `Q_inf` by
[`Q_infinity`](@ref). This function returns a vector with four real
values, the first two are the real and imaginary values of `G` at `ξ₁`
and the second two are their derivatives.
"""
function G(μ::T, γ_real::T, γ_imag::T, κ::T, ϵ::T, ξ₁::T, λ::CGLParams{T}) where {T}
    Q_0, Q_0_dξ = Q_zero(μ, κ, ϵ, ξ₁, λ)

    γ = T == Arb ? Acb(γ_real, γ_imag) : complex(γ_real, γ_imag)
    Q_inf, Q_inf_dξ = Q_infinity(γ, κ, ϵ, ξ₁, λ)

    G1 = Q_0 - Q_inf
    G2 = Q_0_dξ - Q_inf_dξ

    SVector(real(G1), imag(G1), real(G2), imag(G2))
end

"""
    G_jacobian_kappa(μ, γ_real, γ_imag, κ, ϵ, ξ₁, λ::CGLParams)

This function computes the Jacobian of [`G`](@ref) w.r.t. the
parameters `μ`, `γ_real`, `γ_imag` and `κ`.
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
    if λ.d == 3 && iszero(ϵ) && (iswide(μ) || iswide(κ)) && ξ₁ > 100
        μs = mince(μ, ifelse(iswide(μ), 4, 1))
        κs = mince(κ, ifelse(iswide(κ), 96, 1))

        # IMPROVE: We currently fix a tighter tolerance here. This
        # should possibly be adjustable.
        Q_0_Js = tmap(
            ((μ, κ),) -> Q_zero_jacobian_kappa(μ, κ, ϵ, ξ₁, λ, tol = 1e-14),
            collect(Iterators.product(μs, κs)),
        )
        Q_0_J = SMatrix{2,2}(
            Arblib.union(getindex.(Q_0_Js, 1)...),
            Arblib.union(getindex.(Q_0_Js, 2)...),
            Arblib.union(getindex.(Q_0_Js, 3)...),
            Arblib.union(getindex.(Q_0_Js, 4)...),
        )

        γ = T == Arb ? Acb(γ_real, γ_imag) : complex(γ_real, γ_imag)

        Q_inf_Js = tmap(κ -> Q_infinity_jacobian_kappa(γ, κ, ϵ, ξ₁, λ), κs)
        Q_inf_J = SMatrix{2,2}(
            Arblib.union(getindex.(Q_inf_Js, 1)...),
            Arblib.union(getindex.(Q_inf_Js, 2)...),
            Arblib.union(getindex.(Q_inf_Js, 3)...),
            Arblib.union(getindex.(Q_inf_Js, 4)...),
        )
    else
        Q_0_J = Q_zero_jacobian_kappa(μ, κ, ϵ, ξ₁, λ)

        γ = T == Arb ? Acb(γ_real, γ_imag) : complex(γ_real, γ_imag)
        Q_inf_J = Q_infinity_jacobian_kappa(γ, κ, ϵ, ξ₁, λ)
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
    G_jacobian_epsilon(μ, γ_real, γ_imag, κ, ϵ, ξ₁, λ::CGLParams)

This function computes the Jacobian of [`G`](@ref) w.r.t. the
parameters `μ`, `γ_real`, `γ_imag` and `ϵ`.
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
    Q_0_J = Q_zero_jacobian_epsilon(μ, κ, ϵ, ξ₁, λ)

    γ = T == Arb ? Acb(γ_real, γ_imag) : complex(γ_real, γ_imag)
    Q_inf_J = Q_infinity_jacobian_epsilon(γ, κ, ϵ, ξ₁, λ)

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
