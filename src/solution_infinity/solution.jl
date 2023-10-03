"""
    solution_zero(γ, κ, ξ₁, λ::AbstractGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes `[Q(ξ₁), d(Q)(ξ₁)]`.
"""
function solution_infinity(γ::Acb, κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    v = Arb(0.1) # TODO: How to pick this?
    normv = solution_infinity_fixed_point(γ, κ, ξ₁, v, λ)[1]

    # Needed for alternative version of I_P_bound
    normv_dξ =
        C_P_dξ(κ, λ, ξ₁) * abs(γ) * ξ₁^(-v - 1) +
        C_u_dξ(κ, ξ₁, v, λ) * normv^(2λ.σ + 1) * ξ₁^(2λ.σ * v - 1)

    I_E_bound = I_E_0(γ, κ, ξ₁, v, normv, λ)
    I_E_dξ_bound = I_E_dξ_0(γ, κ, ξ₁, v, normv, λ)

    #I_P_bound = I_P_0(γ, κ, ξ₁, v, normv, λ) # Alternative version
    I_P_bound = I_P_0(γ, κ, ξ₁, v, normv, normv_dξ, λ)
    I_P_dξ_bound = I_P_dξ_0(γ, κ, ξ₁, v, normv, λ)

    Q = P(ξ₁, (λ, κ)) * (γ + I_E_bound) + E(ξ₁, (λ, κ)) * I_P_bound

    dQ =
        P(ξ₁, (λ, κ)) * I_E_dξ_bound +
        P_dξ(ξ₁, (λ, κ)) * (γ + I_E_bound) +
        E(ξ₁, (λ, κ)) * I_P_dξ_bound +
        E_dξ(ξ₁, (λ, κ)) * I_P_bound

    return SVector(Q, dQ)
end

function solution_infinity(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::AbstractGLParams{Float64},
)
    # FIXME: Take into account the rest of Q for this case as well?
    Q = γ * P(ξ₁, (λ, κ))
    dQ = γ * P_dξ(ξ₁, (λ, κ))

    return SVector(Q, dQ)
end


"""
    solution_infinity_jacobian(γ, κ, ξ₁, λ::AbstractGLParams)

Let `Q` be the solution to [`fpp_infinity_complex`](@ref). This
function computes Jacobian w.r.t. the parameters `γ` and `κ` of
`[Q(ξ₁), d(Q)(ξ₁)]. The Jacobian is given by
```
[
d(Q(ξ₁), μ) d(Q(ξ₁), κ)
d(d(Q)(ξ₁), μ) d((Q)(ξ₁), κ)
]
```
where we use `d(Q, μ)` to denote the derivative of `Q` w.r.t. `μ`.
"""
function solution_infinity_jacobian(γ::Acb, κ::Arb, ξ₁::Arb, λ::AbstractGLParams{Arb})
    v = Arb(0.1) # TODO: How to pick this?

    normv = solution_infinity_fixed_point(γ, κ, ξ₁, v, λ)[1]

    normv_dξ =
        C_P_dξ(κ, λ, ξ₁) * abs(γ) * ξ₁^(-v - 1) +
        C_u_dξ(κ, ξ₁, v, λ) * normv^(2λ.σ + 1) * ξ₁^(2λ.σ * v - 1)
    normv_dγ = let (CT1, CT2) = C_T1(v, κ, λ, ξ₁), σ = λ.σ
        num = CT1 * ξ₁^-v
        den = (1 - (2σ + 1) * CT2 * ξ₁^(-2 + 2σ * v) * normv^2σ)

        if Arblib.ispositive(den)
            num / den
        else
            @warn "Not positive denominator for normv_dγ" num den
            indeterminate(num)
        end
    end
    normv_dκ = let ξ = ξ₁ # FIXME: We should take supremum for ξ >= ξ₁
        num = (
            C_P_dκ(κ, λ, ξ₁) * abs(γ) * log(ξ) * ξ^-v +
            (C_u_dκ_1(κ, ξ₁, v, λ) * normv + C_u_dκ_2(κ, ξ₁, v, λ) * normv_dξ) * normv^2λ.σ
        )
        den = (1 - C_u_dκ_3(κ, ξ₁, v, λ) * normv^2λ.σ)

        if Arblib.ispositive(den)
            num / den
        else
            @warn "Not positive denominator for normv_dκ" num den
            indeterminate(num)
        end
    end

    I_E_bound = I_E_0(γ, κ, ξ₁, v, normv, λ)
    I_P_bound = I_P_0(γ, κ, ξ₁, v, normv, λ)
    I_E_dξ_bound = I_E_dξ_0(γ, κ, ξ₁, v, normv, λ)
    I_P_dξ_bound = I_P_dξ_0(γ, κ, ξ₁, v, normv, λ)

    I_E_dγ_bound = I_E_dγ_0(γ, κ, ξ₁, v, normv, normv_dγ, λ)
    I_P_dγ_bound = I_P_dγ_0(γ, κ, ξ₁, v, normv, normv_dγ, λ)
    I_E_dξ_dγ_bound = I_E_dξ_dγ_0(γ, κ, ξ₁, v, normv, normv_dγ, λ)
    I_P_dξ_dγ_bound = I_P_dξ_dγ_0(γ, κ, ξ₁, v, normv, normv_dγ, λ)

    I_E_dκ_bound = I_E_dκ_0(γ, κ, ξ₁, v, normv, normv_dκ, λ)
    I_P_dκ_bound = I_P_dκ_0(γ, κ, ξ₁, v, normv, normv_dξ, normv_dκ, λ)
    I_E_dξ_dκ_bound = I_E_dξ_dκ_0(γ, κ, ξ₁, v, normv, normv_dκ, λ)
    I_P_dξ_dκ_bound = I_P_dξ_dκ_0(γ, κ, ξ₁, v, normv, normv_dκ, λ)

    Q_dγ = P(ξ₁, (λ, κ)) * (one(γ) + I_E_dγ_bound) + E(ξ₁, (λ, κ)) * I_P_dγ_bound

    dQ_dγ =
        P(ξ₁, (λ, κ)) * I_E_dξ_dγ_bound +
        P_dξ(ξ₁, (λ, κ)) * (one(γ) + I_E_dγ_bound) +
        E(ξ₁, (λ, κ)) * I_P_dξ_dγ_bound +
        E_dξ(ξ₁, (λ, κ)) * I_P_dγ_bound

    Q_dκ =
        P(ξ₁, (λ, κ)) * I_E_dκ_bound +
        P_dκ(ξ₁, (λ, κ)) * (γ + I_E_bound) +
        E(ξ₁, (λ, κ)) * I_P_dκ_bound +
        E_dκ(ξ₁, (λ, κ)) * I_P_bound

    dQ_dκ =
        P(ξ₁, (λ, κ)) * I_E_dξ_dκ_bound +
        P_dξ(ξ₁, (λ, κ)) * I_E_dκ_bound +
        P_dκ(ξ₁, (λ, κ)) * I_E_dξ_bound +
        P_dξ_dκ(ξ₁, (λ, κ)) * (γ + I_E_bound) +
        E(ξ₁, (λ, κ)) * I_P_dξ_dκ_bound +
        E_dξ(ξ₁, (λ, κ)) * I_P_dκ_bound +
        E_dκ(ξ₁, (λ, κ)) * I_P_dξ_bound +
        E_dξ_dκ(ξ₁, (λ, κ)) * I_P_bound

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end

function solution_infinity_jacobian(
    γ::ComplexF64,
    κ::Float64,
    ξ₁::Float64,
    λ::AbstractGLParams{Float64},
)
    # FIXME: Take into account the rest of Q
    Q_dγ = P(ξ₁, (λ, κ))
    dQ_dγ = P_dξ(ξ₁, (λ, κ))

    Q_dκ = γ * P_dκ(ξ₁, (λ, κ))
    dQ_dκ = γ * P_dξ_dκ(ξ₁, (λ, κ))

    return SMatrix{2,2}(Q_dγ, dQ_dγ, Q_dκ, dQ_dκ)
end
