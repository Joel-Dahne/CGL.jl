"""
    I_E_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, norm_u::Arb, λ::AbstractGLParams{Arb})

Let
```
I_E(ξ) = ∫_ξ₁^ξ J_E(η) * abs(u(η))^2σ * u(η) dη
```
This computes a complex ball enclosing `I_E(ξ₁)`.

Note that that `I_E(ξ₁)` is exactly zero. We keep the method mostly
for consistency with the other methods.
"""
function I_E_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, norm_u::Arb, λ::AbstractGLParams{Arb})
    return zero(γ)
end

"""
    I_P_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, norm_u::Arb, λ::AbstractGLParams{Arb})

Let
```
I_P(ξ) = ∫_ξ^∞ J_P(η) * abs(u(η))^2σ * u(η) dη
```
This computes a complex ball enclosing `I_P(ξ₁)`.
"""
function I_P_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, norm_u::Arb, λ::AbstractGLParams{Arb})
    bound = I_P_0_bound_1(κ, ξ₁, v, norm_u, λ)

    return add_error(zero(γ), bound)
end

function I_P_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    λ::AbstractGLParams{Arb},
)
    # The separate bounds are given in the paper. We compute both and
    # take the minimum.
    bound = min(
        I_P_0_bound_1(κ, ξ₁, v, norm_u, λ),
        I_P_0_bound_2(κ, ξ₁, v, norm_u, norm_u_dξ, λ),
    )

    return add_error(zero(γ), bound)
end

function I_P_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    λ::AbstractGLParams{Arb},
)
    # The separate bounds are given in the paper. We compute both and
    # take the minimum.
    bound = min(
        I_P_0_bound_1(κ, ξ₁, v, norm_u, λ),
        I_P_0_bound_2(κ, ξ₁, v, norm_u, norm_u_dξ, λ),
        I_P_0_bound_3(κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, λ),
    )

    return add_error(zero(γ), bound)
end

function I_P_0_bound_1(κ::Arb, ξ₁::Arb, v::Arb, norm_u::Arb, λ::AbstractGLParams{Arb})
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0 # Required for integral to converge

    bound =
        C_I_P(κ, ξ₁, v, λ) *
        norm_u^(2σ + 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    return bound
end

# Similar to the above method but doing one step of partial integration
function I_P_0_bound_2(
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < 0 # Required for integral to converge

    bound =
        (C_I_P_1_1(κ, ξ₁, v, λ) * norm_u * ξ₁^(-1) + C_I_P_1_2(κ, ξ₁, v, λ) * norm_u_dξ) *
        norm_u^2σ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    return bound
end

# Similar to the above method but doing two steps of partial integration
function I_P_0_bound_3(
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 4 < 0 # Required for integral to converge

    bound =
        (
            C_I_P_2_1(κ, ξ₁, v, λ) * norm_u^2 +
            C_I_P_2_2(κ, ξ₁, v, λ) * norm_u^2 * ξ₁^-2 +
            C_I_P_2_3(κ, ξ₁, v, λ) * norm_u * norm_u_dξ * ξ₁^-1 +
            C_I_P_2_4(κ, ξ₁, v, λ) * norm_u_dξ^2 +
            C_I_P_2_5(κ, ξ₁, v, λ) * norm_u * norm_u_dξ_dξ
        ) *
        norm_u^(2σ - 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 4)

    return bound
end

function I_E_dξ_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, norm_u::Arb, λ::AbstractGLParams{Arb})
    (; σ) = λ

    bound = C_I_E_dξ(κ, ξ₁, v, λ) * norm_u^(2σ + 1) * ξ₁^((2σ + 1) * v - 3)

    return add_error(zero(γ), bound)
end

function I_P_dξ_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, norm_u::Arb, λ::AbstractGLParams{Arb})
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < -1 # Required for integral to converge

    bound =
        C_I_P_dξ(κ, ξ₁, v, λ) *
        norm_u^(2σ + 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    return add_error(zero(γ), bound)
end

###
# Derivatives with respect to γ
###

function I_E_dγ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    return zero(γ) # This is exactly zero
end

function I_P_dγ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    bound = I_P_dγ_0_bound_1(κ, ξ₁, v, norm_u, norm_u_dγ, λ)

    return add_error(zero(γ), bound)
end

function I_P_dγ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    # The separate bounds are given in the paper. We compute both and
    # take the minimum.
    bound = min(
        I_P_dγ_0_bound_1(κ, ξ₁, v, norm_u, norm_u_dγ, λ),
        I_P_dγ_0_bound_2(κ, ξ₁, v, norm_u, norm_u_dγ, norm_u_dξ, norm_u_dξ_dγ, λ),
    )

    return add_error(zero(γ), bound)
end

function I_P_dγ_0_bound_1(
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0 # Required for integral to converge

    bound =
        (2σ + 1) *
        C_I_P(κ, ξ₁, v, λ) *
        norm_u^2σ *
        norm_u_dγ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    return bound
end

# Similar to the above method but doing one step of partial integration
function I_P_dγ_0_bound_2(
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < 0 # Required for integral to converge

    bound =
        (
            (2σ + 1) * C_I_P_1_1(κ, ξ₁, v, λ) * norm_u * norm_u_dγ * ξ₁^(-1) +
            C_I_P_1_2(κ, ξ₁, v, λ) * (2σ * norm_u_dξ * norm_u_dγ + norm_u_dξ_dγ)
        ) *
        norm_u^(2σ - 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    return bound
end
function I_E_dξ_dγ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; σ) = λ

    bound = (2σ + 1) * C_I_E_dξ(κ, ξ₁, v, λ) * norm_u^2σ * norm_u_dγ * ξ₁^((2σ + 1) * v - 3)

    return add_error(zero(γ), bound)
end

function I_P_dξ_dγ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < -1 # Required for integral to converge

    bound =
        (2σ + 1) *
        C_I_P_dξ(κ, ξ₁, v, λ) *
        norm_u^2σ *
        norm_u_dγ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    return add_error(zero(γ), bound)
end

function I_E_dκ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    return zero(γ) # This is exactly zero
end

function I_P_dκ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    norm_u_dκ::Arb,
    norm_u_dξ_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    return I_P_dκ_1_0(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dξ_dξ, norm_u_dκ, λ) +
           I_P_dκ_2_0(γ, κ, ξ₁, v, norm_u, norm_u_dξ, norm_u_dκ, norm_u_dξ_dκ, λ)
end

function I_P_dκ_1_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dξ_dξ::Arb,
    norm_u_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    bound1 =
        C_I_P_dκ_1(κ, ξ₁, v, λ) *
        norm_u^(2σ + 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    bound2 =
        (C_I_P_dκ_2(κ, ξ₁, v, λ) * norm_u * ξ₁^-1 + C_I_P_dκ_4(κ, ξ₁, v, λ) * norm_u_dξ) *
        norm_u^2σ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    bound3 =
        (
            C_I_P_dκ_3(κ, ξ₁, v, λ) * norm_u^2 * ξ₁^-2 +
            C_I_P_dκ_5(κ, ξ₁, v, λ) * norm_u * norm_u_dξ * ξ₁^-1 +
            C_I_P_dκ_6(κ, ξ₁, v, λ) * norm_u_dξ^2 +
            C_I_P_dκ_7(κ, ξ₁, v, λ) * norm_u * norm_u_dξ_dξ
        ) *
        norm_u^(2σ - 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    bound = bound1 + bound2 + bound3

    return add_error(zero(γ), bound)
end

function I_P_dκ_2_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dξ::Arb,
    norm_u_dκ::Arb,
    norm_u_dξ_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    bound =
        (
            (2σ + 1) * C_I_P_1_1(κ, ξ₁, v, λ) * norm_u * norm_u_dκ * ξ₁^(-1) +
            C_I_P_1_2(κ, ξ₁, v, λ) * (2σ * norm_u_dξ * norm_u_dκ + norm_u_dξ_dκ)
        ) *
        norm_u^(2σ - 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    return add_error(zero(γ), bound)
end

function I_E_dξ_dκ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; σ) = λ

    bound =
        (
            C_I_E_dξ_dκ_1(κ, ξ₁, v, λ) * norm_u * log(ξ₁) +
            C_I_E_dξ_dκ_2(κ, ξ₁, v, λ) * norm_u_dκ
        ) *
        norm_u^2σ *
        ξ₁^((2σ + 1) * v - 3)

    return add_error(zero(γ), bound)
end

function I_P_dξ_dκ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    norm_u::Arb,
    norm_u_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ, ϵ) = λ

    _, _, c = _abc(κ, λ)

    bound =
        (
            C_I_P_dξ_dκ_1(κ, ξ₁, v, λ) * norm_u +
            C_I_P_dξ_dκ_2(κ, ξ₁, v, λ) * norm_u_dκ * ξ₁^(-2)
        ) *
        norm_u^2σ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 1)

    return add_error(zero(γ), bound)
end
