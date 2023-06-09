"""
    enclose_derivative_F(κ::Arb, μ::Arb, p::AbstractGLParams{Arb}, ξ₁::Arb, v::Arb; non_rigorous::Bool = false)

Compute an enclosure of the derivative with respect to `ξ` of
```
u_0(ξ, (p, κ, μ)) - u_T(ξ, (p, κ, γ))
```
evaluated at `ξ = ξ₁`.

Here `u_0(ξ, (p, κ, μ))` is the solution satisfying the initial
condition `u_0(0, (p, κ, μ)) = μ` and `u_T(ξ, (p, κ, γ))` is the fix
point of the operator `T` with `γ` taken as in
[`check_existence_fixed_point`](@ref).
"""
function enclose_derivative_F(
    κ::Arb,
    μ::Arb,
    p::AbstractGLParams{Arb},
    ξ₁::Arb,
    v::Arb;
    non_rigorous::Bool = false,
)
    # Compute u_0(ξ₁, (p, κ, μ))
    u_0, u_0_dξ = let
        sol = if non_rigorous
            solution_zero_float(κ, μ, p, ξ₁)
        else
            solution_zero_capd(κ, μ, p, ξ₁)
        end
        Acb(sol[1], sol[2]), Acb(sol[3], sol[4])
    end

    if non_rigorous
        # Add some extra radius to see if it still works
        u_0 = add_error(u_0, max(radius(real(u_0)), radius(imag(u_0))))
        u_0_dξ = add_error(u_0_dξ, max(radius(real(u_0_dξ)), radius(imag(u_0_dξ))))
    end

    # Compute enclosure of γ
    γ = solution_infinity_γ(κ, μ, p, ξ₁, v, u_0)

    u_T_dξ = let
        C₁, C₂ = C_fix_point(abs(γ), v, κ, p, ξ₁)
        C₃, C₄ = C_fix_point_dξ(abs(γ), v, κ, p, ξ₁)

        real_c = κ * p.ϵ / (1 + p.ϵ)^2

        f_bound =
            add_error(zero(C₁), C₁ * (ξ₁^(2p.σ * v + v - 2) - ξ₁^(2p.σ * v + v - 2)))
        g_bound = add_error(
            zero(C₂),
            C₂ * exp(-real_c * ξ₁^2) * ξ₁^(-2 / p.σ + 2p.σ * v + v + p.d - 2),
        )

        f_dξ_bound = add_error(zero(C₃), C₃ * ξ₁^(2p.σ * v + v - 3))
        g_dξ_bound = add_error(
            zero(C₄),
            C₄ * exp(-real_c * ξ₁^2) * ξ₁^(-2 / p.σ + 2p.σ * v + v + p.d - 3),
        )

        # IMPROVE: Optimize these enclosures in κ
        P_ξ₁, P_dξ_ξ₁ = let tmp = P(ArbSeries((ξ₁, 1)), (p, κ))
            tmp[0], tmp[1]
        end
        E_ξ₁, E_dξ_ξ₁ = let tmp = E(ArbSeries((ξ₁, 1)), (p, κ))
            tmp[0], tmp[1]
        end

        γ * P_dξ_ξ₁ +
        P_dξ_ξ₁ * f_bound +
        P_ξ₁ * f_dξ_bound +
        E_dξ_ξ₁ * g_bound +
        E_ξ₁ * g_dξ_bound
    end

    return u_0_dξ - u_T_dξ
end
