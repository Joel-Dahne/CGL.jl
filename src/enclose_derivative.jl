"""
    enclose_derivative_F(κ::Arb, μ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb, v::Arb; non_rigorous::Bool = false)

Let `Q_0` be the solution to [`ivp_zero_complex`](@ref) and `Q_inf` a
solution to [`fpp_infinity_complex`](@ref). This function computes an
enclosure of the derivative with respect to `ξ` of
```
d(Q_0)(ξ₁) - d(Q_inf)(ξ₁)
```
We here use `d(Q)` to denote the derivative of `Q` with respect to
`ξ`.
"""
function enclose_derivative_F(
    μ::Arb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::AbstractGLParams{Arb};
    non_rigorous::Bool = false,
)
    # Compute Q_0(ξ₁)
    Q_0, Q_0_dξ = if non_rigorous
        let sol = solution_zero_float(μ, κ, ξ₁, λ)
            Acb(sol[1], sol[2]), Acb(sol[3], sol[4])
        end
    else
        solution_zero(μ, κ, ξ₁, λ)
    end

    if non_rigorous
        # Add some extra radius to see if it still works
        Q_0 = add_error(Q_0, max(radius(real(Q_0)), radius(imag(Q_0))))
        Q_0_dξ = add_error(Q_0_dξ, max(radius(real(Q_0_dξ)), radius(imag(Q_0_dξ))))
    end

    # Compute enclosure of γ
    γ = solution_infinity_γ(κ, μ, λ, ξ₁, v, Q_0)

    Q_inf_dξ = let
        C₁, C₂ = C_fix_point(abs(γ), v, κ, λ, ξ₁)
        C₃, C₄ = C_fix_point_dξ(abs(γ), v, κ, λ, ξ₁)

        real_c = κ * λ.ϵ / (1 + λ.ϵ)^2

        f_bound =
            add_error(zero(C₁), C₁ * (ξ₁^(2λ.σ * v + v - 2) - ξ₁^(2λ.σ * v + v - 2)))
        g_bound = add_error(
            zero(C₂),
            C₂ * exp(-real_c * ξ₁^2) * ξ₁^(-2 / λ.σ + 2λ.σ * v + v + λ.d - 2),
        )

        f_dξ_bound = add_error(zero(C₃), C₃ * ξ₁^(2λ.σ * v + v - 3))
        g_dξ_bound = add_error(
            zero(C₄),
            C₄ * exp(-real_c * ξ₁^2) * ξ₁^(-2 / λ.σ + 2λ.σ * v + v + λ.d - 3),
        )

        # IMPROVE: Optimize these enclosures in κ
        P_ξ₁, P_dξ_ξ₁ = let tmp = P(ArbSeries((ξ₁, 1)), (λ, κ))
            tmp[0], tmp[1]
        end
        E_ξ₁, E_dξ_ξ₁ = let tmp = E(ArbSeries((ξ₁, 1)), (λ, κ))
            tmp[0], tmp[1]
        end

        γ * P_dξ_ξ₁ +
        P_dξ_ξ₁ * f_bound +
        P_ξ₁ * f_dξ_bound +
        E_dξ_ξ₁ * g_bound +
        E_ξ₁ * g_dξ_bound
    end

    return Q_0_dξ - Q_inf_dξ
end
