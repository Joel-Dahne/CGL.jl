struct NormBounds
    Q::Arb
    Q_dξ::Arb
    Q_dξ_dξ::Arb
    Q_dξ_dξ_dξ::Arb
    Q_dγ::Arb
    Q_dγ_dξ::Arb
    Q_dκ::Arb
    Q_dκ_dξ::Arb
    Q_dϵ::Arb
    Q_dϵ_dξ::Arb
    Q2σQ::Arb
    Q2σQ_dξ::Arb
    Q2σQ_dξ_dξ::Arb
    Q2σQ_dξ_dξ_dξ::Arb
    Q2σQ_dγ::Arb
    Q2σQ_dγ_dξ::Arb
    Q2σQ_dκ::Arb
    Q2σQ_dκ_dξ::Arb
    Q2σQ_dϵ::Arb
    Q2σQ_dϵ_dξ::Arb

    function NormBounds(
        γ::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        v::Arb,
        λ::CGLParams{Arb},
        C::FunctionBounds;
        include_dκ::Bool = false,
        include_dϵ::Bool = false,
    )
        (; σ) = λ

        norms = new(
            norm_bound_Q(γ, κ, ϵ, ξ₁, v, λ, C),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
        )

        norms.Q_dξ[] = norm_bound_Q_dξ(γ, κ, ϵ, ξ₁, v, λ, C, norms)
        norms.Q_dξ_dξ[] = norm_bound_Q_dξ_dξ(γ, κ, ϵ, ξ₁, v, λ, C, norms)
        norms.Q_dξ_dξ_dξ[] = norm_bound_Q_dξ_dξ_dξ(γ, κ, ϵ, ξ₁, v, λ, C, norms)

        # Norms of abs(Q)^2σ * Q and its derivatives
        norm_Q_series =
            ArbSeries((norms.Q, norms.Q_dξ, norms.Q_dξ_dξ / 2, norms.Q_dξ_dξ_dξ / 6))
        norm_Q2σQ_series = norm_Q_series^2σ * norm_Q_series

        norms.Q2σQ[] = norm_Q2σQ_series[0]
        norms.Q2σQ_dξ[] = norm_Q2σQ_series[1]
        norms.Q2σQ_dξ_dξ[] = 2norm_Q2σQ_series[2]
        norms.Q2σQ_dξ_dξ_dξ[] = 6norm_Q2σQ_series[3]

        if include_dκ || include_dϵ
            norms.Q_dγ[] = norm_bound_Q_dγ(γ, κ, ϵ, ξ₁, v, λ, C, norms)
            norms.Q_dγ_dξ[] = norm_bound_Q_dγ_dξ(γ, κ, ϵ, ξ₁, v, λ, C, norms)

            # Norms of abs(Q)^2σ * Q differentiated w.r.t. γ
            norms.Q2σQ_dγ[] = (2σ + 1) * norms.Q^2σ * norms.Q_dγ
            norms.Q2σQ_dγ_dξ[] =
                (2σ + 1) *
                norms.Q^(2σ - 1) *
                (2σ * norms.Q_dξ * norms.Q_dγ + norms.Q * norms.Q_dγ_dξ)
        end

        if include_dκ
            norms.Q_dκ[] = norm_bound_Q_dκ(γ, κ, ϵ, ξ₁, v, λ, C, norms)
            norms.Q_dκ_dξ[] = norm_bound_Q_dκ_dξ(γ, κ, ϵ, ξ₁, v, λ, C, norms)

            # Norms of abs(Q)^2σ * Q differentiated w.r.t. κ
            norms.Q2σQ_dκ[] = (2σ + 1) * norms.Q^2σ * norms.Q_dκ
            norms.Q2σQ_dκ_dξ[] =
                (2σ + 1) *
                norms.Q^(2σ - 1) *
                (2σ * norms.Q_dξ * norms.Q_dκ + norms.Q * norms.Q_dκ_dξ)
        end

        if include_dϵ
            norms.Q_dϵ[] = norm_bound_Q_dϵ(γ, κ, ϵ, ξ₁, v, λ, C, norms)
            norms.Q_dϵ_dξ[] = norm_bound_Q_dϵ_dξ(γ, κ, ϵ, ξ₁, v, λ, C, norms)

            # Norms of abs(Q)^2σ * Q differentiated w.r.t. ϵ
            norms.Q2σQ_dϵ[] = (2σ + 1) * norms.Q^2σ * norms.Q_dϵ
            norms.Q2σQ_dϵ_dξ[] =
                (2σ + 1) *
                norms.Q^(2σ - 1) *
                (2σ * norms.Q_dξ * norms.Q_dϵ + norms.Q * norms.Q_dϵ_dξ)
        end

        return norms
    end
end

norm_bound_Q(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
) = solution_infinity_fixed_point(γ, κ, ϵ, ξ₁, v, λ, C)[1]

function norm_bound_Q_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    norms::NormBounds,
)
    (; σ) = λ
    return C.P_dξ * abs(γ) * ξ₁^(-v - 1) +
           C_u_dξ(κ, ϵ, ξ₁, v, λ, C) * norms.Q^(2σ + 1) * ξ₁^(2σ * v - 1)
end

function norm_bound_Q_dξ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    norms::NormBounds,
)
    (; σ) = λ
    return C.P_dξ_dξ * abs(γ) * ξ₁^(-v - 2) +
           (
               C_u_dξ_dξ_1(κ, ϵ, ξ₁, v, λ, C) * norms.Q * ξ₁^(-1) +
               C_u_dξ_dξ_2(κ, ϵ, ξ₁, v, λ, C) * norms.Q_dξ
           ) *
           norms.Q^2σ *
           ξ₁^(2σ * v - 1)
end

# TODO: It seems like norm_Q_dξ_dξ_dξ is larger than norm_Q_dξ_dξ.
# Previously the norms have been decreasing for higher
# derivatives. See if this makes sense.
function norm_bound_Q_dξ_dξ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    norms::NormBounds,
)
    (; σ) = λ
    return C.P_dξ_dξ_dξ * abs(γ) * ξ₁^(-v - 3) +
           (
               C_u_dξ_dξ_dξ_1(κ, ϵ, ξ₁, v, λ, C) * norms.Q^2 +
               C_u_dξ_dξ_dξ_2(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ * ξ₁^(-1) +
               C_u_dξ_dξ_dξ_3(κ, ϵ, ξ₁, v, λ, C) * norms.Q_dξ^2 +
               C_u_dξ_dξ_dξ_4(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ_dξ
           ) *
           norms.Q^(2σ - 1) *
           ξ₁^(2σ * v - 1)
end

function norm_bound_Q_dγ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    norms::NormBounds,
)
    (; σ) = λ
    num = C.P * ξ₁^-v
    den = (1 - (2σ + 1) * C_T1(κ, ϵ, ξ₁, v, λ, C) * ξ₁^(-2 + 2σ * v) * norms.Q^2σ)

    if Arblib.ispositive(den)
        num / den
    else
        indeterminate(num)
    end
end

function norm_bound_Q_dγ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    norms::NormBounds,
)
    (; σ) = λ
    return C.P_dξ * ξ₁^(-v - 1) +
           (2σ + 1) *
           C_u_dξ(κ, ϵ, ξ₁, v, λ, C) *
           norms.Q^2σ *
           norms.Q_dγ *
           ξ₁^(2λ.σ * v - 1)
end

function norm_bound_Q_dκ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    norms::NormBounds,
)
    (; σ) = λ
    num = (
        C_u_dκ_1(κ, ϵ, ξ₁, v, λ, C) * abs(γ) +
        (
            C_u_dκ_2(κ, ϵ, ξ₁, v, λ, C) * norms.Q^2 +
            C_u_dκ_3(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ +
            C_u_dκ_4(κ, ϵ, ξ₁, v, λ, C) * norms.Q_dξ^2 +
            C_u_dκ_5(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ_dξ
        ) * norms.Q^(2σ - 1)
    )
    den = (1 - C_u_dκ_6(κ, ϵ, ξ₁, v, λ, C) * norms.Q^2σ)

    if Arblib.ispositive(den)
        num / den
    else
        indeterminate(num)
    end
end

function norm_bound_Q_dκ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    norms::NormBounds,
)
    (; σ) = λ
    @assert ξ₁ >= ℯ
    return C.P_dκ * abs(γ) * log(ξ₁) * ξ₁^(-v - 1) +
           (
        C_u_dξ_dκ_1(κ, ϵ, ξ₁, v, λ, C) * norms.Q^2 +
        C_u_dξ_dκ_2(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dκ +
        C_u_dξ_dκ_3(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ +
        C_u_dξ_dκ_4(κ, ϵ, ξ₁, v, λ, C) * norms.Q_dξ^2 +
        C_u_dξ_dκ_5(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ_dξ
    ) * norms.Q^(2σ - 1)
end

function norm_bound_Q_dϵ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    norms::NormBounds,
)
    (; σ) = λ
    num = (
        C_u_dϵ_1(κ, ϵ, ξ₁, v, λ, C) * abs(γ) +
        (
            C_u_dϵ_2(κ, ϵ, ξ₁, v, λ, C) * norms.Q^2 +
            C_u_dϵ_3(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ +
            C_u_dϵ_4(κ, ϵ, ξ₁, v, λ, C) * norms.Q_dξ^2 +
            C_u_dϵ_5(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ_dξ
        ) * norms.Q^(2σ - 1)
    )
    den = (1 - C_u_dϵ_6(κ, ϵ, ξ₁, v, λ, C) * norms.Q^2σ)

    if Arblib.ispositive(den)
        num / den
    else
        indeterminate(num)
    end
end

function norm_bound_Q_dϵ_dξ(
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds,
    norms::NormBounds,
)
    (; σ) = λ

    return C.P_dϵ * abs(γ) * ξ₁^(-v - 1) +
           (
        C_u_dξ_dϵ_1(κ, ϵ, ξ₁, v, λ, C) * norms.Q^2 +
        C_u_dξ_dϵ_2(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dϵ +
        C_u_dξ_dϵ_3(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ +
        C_u_dξ_dϵ_4(κ, ϵ, ξ₁, v, λ, C) * norms.Q_dξ^2 +
        C_u_dξ_dϵ_5(κ, ϵ, ξ₁, v, λ, C) * norms.Q * norms.Q_dξ_dξ
    ) * norms.Q^(2σ - 1)
end
