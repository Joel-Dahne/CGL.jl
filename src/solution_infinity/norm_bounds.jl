norm_bound_u(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
) = solution_infinity_fixed_point(γ, κ, ϵ, ξ₁, v, λ, C)[1]

function norm_bound_u_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; σ) = λ
    return C.P_dξ * abs(γ) * ξ₁^(-v - 1) +
           C_u_dξ(κ, ϵ, ξ₁, v, λ, C) * norm_u^(2σ + 1) * ξ₁^(2σ * v - 1)
end

function norm_bound_u_dξ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; σ) = λ
    return C.P_dξ_dξ * abs(γ) * ξ₁^(-v - 2) +
           (
               C_u_dξ_dξ_1(κ, ϵ, ξ₁, v, λ, C) * norm_u * ξ₁^(-1) +
               C_u_dξ_dξ_2(κ, ϵ, ξ₁, v, λ, C) * norm_u_dξ
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
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; σ) = λ
    return C.P_dξ_dξ_dξ * abs(γ) * ξ₁^(-v - 3) +
           (
               C_u_dξ_dξ_dξ_1(κ, ϵ, ξ₁, v, λ, C) * norm_u^2 +
               C_u_dξ_dξ_dξ_2(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dξ * ξ₁^(-1) +
               C_u_dξ_dξ_dξ_3(κ, ϵ, ξ₁, v, λ, C) * norm_u_dξ^2 +
               C_u_dξ_dξ_dξ_4(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dξ_dξ
           ) *
           norm_u^(2σ - 1) *
           ξ₁^(2σ * v - 1)
end

function norm_bound_u_dγ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; σ) = λ
    num = C.P * ξ₁^-v
    den = (1 - (2σ + 1) * C_T1(κ, ϵ, ξ₁, v, λ, C) * ξ₁^(-2 + 2σ * v) * norm_u^2σ)

    if Arblib.ispositive(den)
        num / den
    else
        indeterminate(num)
    end
end

function norm_bound_u_dξ_dγ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; σ) = λ
    return C.P_dξ * ξ₁^(-v - 1) +
           (2σ + 1) * C_u_dξ(κ, ϵ, ξ₁, v, λ, C) * norm_u^2σ * norm_u_dγ * ξ₁^(2λ.σ * v - 1)
end

function norm_bound_u_dκ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; σ) = λ
    num = (
        C_u_dκ_1(κ, ϵ, ξ₁, v, λ, C) * abs(γ) +
        (
            C_u_dκ_2(κ, ϵ, ξ₁, v, λ, C) * norm_u^2 +
            C_u_dκ_3(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dξ +
            C_u_dκ_4(κ, ϵ, ξ₁, v, λ, C) * norm_u_dξ^2 +
            C_u_dκ_5(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dξ_dξ
        ) * norm_u^(2σ - 1)
    )
    den = (1 - C_u_dκ_6(κ, ϵ, ξ₁, v, λ, C) * norm_u^2σ)

    if Arblib.ispositive(den)
        num / den
    else
        indeterminate(num)
    end
end

function norm_bound_u_dξ_dκ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    norm_u_dκ::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; σ) = λ
    @assert ξ₁ >= ℯ
    return C.P_dκ * abs(γ) * log(ξ₁) * ξ₁^(-v - 1) +
           (
        C_u_dξ_dκ_1(κ, ϵ, ξ₁, v, λ, C) * norm_u^2 +
        C_u_dξ_dκ_2(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dκ +
        C_u_dξ_dκ_3(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dξ +
        C_u_dξ_dκ_4(κ, ϵ, ξ₁, v, λ, C) * norm_u_dξ^2 +
        C_u_dξ_dκ_5(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dξ_dξ
    ) * norm_u^(2σ - 1)
end

function norm_bound_u_dϵ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; σ) = λ
    num = (
        C_u_dϵ_1(κ, ϵ, ξ₁, v, λ, C) * abs(γ) +
        (
            C_u_dϵ_2(κ, ϵ, ξ₁, v, λ, C) * norm_u^2 +
            C_u_dϵ_3(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dξ +
            C_u_dϵ_4(κ, ϵ, ξ₁, v, λ, C) * norm_u_dξ^2 +
            C_u_dϵ_5(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dξ_dξ
        ) * norm_u^(2σ - 1)
    )
    den = (1 - C_u_dϵ_6(κ, ϵ, ξ₁, v, λ, C) * norm_u^2σ)

    if Arblib.ispositive(den)
        num / den
    else
        indeterminate(num)
    end
end

function norm_bound_u_dξ_dϵ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    norm_u_dϵ::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
)
    (; σ) = λ

    return C.P_dϵ * abs(γ) * ξ₁^(-v - 1) +
           (
        C_u_dξ_dϵ_1(κ, ϵ, ξ₁, v, λ, C) * norm_u^2 +
        C_u_dξ_dϵ_2(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dϵ +
        C_u_dξ_dϵ_3(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dξ +
        C_u_dξ_dϵ_4(κ, ϵ, ξ₁, v, λ, C) * norm_u_dξ^2 +
        C_u_dξ_dϵ_5(κ, ϵ, ξ₁, v, λ, C) * norm_u * norm_u_dξ_dξ
    ) * norm_u^(2σ - 1)
end
