norm_bound_u(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb}) =
    solution_infinity_fixed_point(γ, κ, ξ₁, v, λ)[1]

function norm_bound_u_dξ(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, norm_u::Arb, λ::CGLParams{Arb})
    (; σ) = λ
    return C_P_dξ(κ, λ, ξ₁) * abs(γ) * ξ₁^(-v - 1) +
           C_u_dξ(κ, ξ₁, v, λ) * norm_u^(2σ + 1) * ξ₁^(2σ * v - 1)
end

function norm_bound_u_dξ_dξ(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    λ::CGLParams{Arb},
)
    (; σ) = λ
    return C_P_dξ_dξ(κ, λ, ξ₁) * abs(γ) * ξ₁^(-v - 2) +
           (
               C_u_dξ_dξ_1(κ, ξ₁, v, λ) * norm_u * ξ₁^(-1) +
               C_u_dξ_dξ_2(κ, ξ₁, v, λ) * norm_u_dξ
           ) *
           norm_u^2σ *
           ξ₁^(2σ * v - 1)
end

# TODO: It seems like norm_u_dξ_dξ_dξ is larger than norm_u_dξ_dξ.
# Previously the norms have been decreasing for higher
# derivatives. See if this makes sense.
function norm_bound_u_dξ_dξ_dξ(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    λ::CGLParams{Arb},
)
    (; σ) = λ
    return C_P_dξ_dξ_dξ(κ, λ, ξ₁) * abs(γ) * ξ₁^(-v - 3) +
           (
               C_u_dξ_dξ_dξ_1(κ, ξ₁, v, λ) * norm_u^2 +
               C_u_dξ_dξ_dξ_2(κ, ξ₁, v, λ) * norm_u * norm_u_dξ * ξ₁^(-1) +
               C_u_dξ_dξ_dξ_3(κ, ξ₁, v, λ) * norm_u_dξ^2 +
               C_u_dξ_dξ_dξ_4(κ, ξ₁, v, λ) * norm_u * norm_u_dξ_dξ
           ) *
           norm_u^(2σ - 1) *
           ξ₁^(2σ * v - 1)
end

function norm_bound_u_dγ(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, norm_u::Arb, λ::CGLParams{Arb})
    (; σ) = λ
    CT1, CT2 = C_T1(v, κ, λ, ξ₁)
    num = CT1 * ξ₁^-v
    den = (1 - (2σ + 1) * CT2 * ξ₁^(-2 + 2σ * v) * norm_u^2σ)

    if Arblib.ispositive(den)
        num / den
    else
        @warn "Not positive denominator when solving for norm_u_dγ" num den
        indeterminate(num)
    end
end

function norm_bound_u_dξ_dγ(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    λ::CGLParams{Arb},
)
    (; σ) = λ
    return C_P_dξ(κ, λ, ξ₁) * ξ₁^(-v - 1) +
           (2σ + 1) * C_u_dξ(κ, ξ₁, v, λ) * norm_u^2σ * norm_u_dγ * ξ₁^(2λ.σ * v - 1)
end

function norm_bound_u_dκ(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    λ::CGLParams{Arb},
)
    (; σ) = λ
    num = (
        C_u_dκ_1(κ, ξ₁, v, λ) * abs(γ) +
        (
            C_u_dκ_2(κ, ξ₁, v, λ) * norm_u^2 +
            C_u_dκ_3(κ, ξ₁, v, λ) * norm_u * norm_u_dξ +
            C_u_dκ_4(κ, ξ₁, v, λ) * norm_u_dξ^2 +
            C_u_dκ_5(κ, ξ₁, v, λ) * norm_u * norm_u_dξ_dξ
        ) * norm_u^(2σ - 1)
    )
    den = (1 - C_u_dκ_6(κ, ξ₁, v, λ) * norm_u^2σ)

    if Arblib.ispositive(den)
        num / den
    else
        @warn "Not positive denominator when solving for norm_u_dγ" num den
        indeterminate(num)
    end
end

function norm_bound_u_dξ_dκ(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    norm_u_dκ::Arb,
    λ::CGLParams{Arb},
)
    (; σ) = λ
    return C_P_dκ(κ, λ, ξ₁) * abs(γ) * log(ξ₁) * ξ₁^(-v - 1) +
           (
        C_u_dξ_dκ_1(κ, ξ₁, v, λ) * norm_u^2 +
        C_u_dξ_dκ_2(κ, ξ₁, v, λ) * norm_u * norm_u_dκ +
        C_u_dξ_dκ_3(κ, ξ₁, v, λ) * norm_u * norm_u_dξ +
        C_u_dξ_dκ_4(κ, ξ₁, v, λ) * norm_u_dξ^2 +
        C_u_dξ_dκ_5(κ, ξ₁, v, λ) * norm_u * norm_u_dξ_dξ
    ) * norm_u^(2σ - 1)
end
