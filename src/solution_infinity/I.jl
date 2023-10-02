"""
    I_E_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, normv::Arb, λ::AbstractGLParams{Arb})

Let
```
I_E(ξ) = ∫_ξ₁^ξ J_E(η) * abs(u(η))^2σ * u(η) dη
```
This computes a complex ball enclosing `I_E(ξ₁)`.

Note that that `I_E(ξ₁)` is exactly zero. We keep the method mostly
for consistency with the other methods.
"""
function I_E_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, normv::Arb, λ::AbstractGLParams{Arb})
    return zero(γ)
end

"""
    I_P_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, normv::Arb, λ::AbstractGLParams{Arb})

Let
```
I_P(ξ) = ∫_ξ^∞ J_P(η) * abs(u(η))^2σ * u(η) dη
```
This computes a complex ball enclosing `I_P(ξ₁)`.
"""
function I_P_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, normv::Arb, λ::AbstractGLParams{Arb})
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0 # Required for integral to converge

    bound =
        C_I_P(κ, ξ₁, v, λ) *
        normv^(2σ + 1) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    return add_error(zero(γ), bound)
end

function I_E_dξ_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, normv::Arb, λ::AbstractGLParams{Arb})
    (; σ) = λ

    bound = C_I_E_dξ(κ, ξ₁, v, λ) * normv^(2σ + 1) * ξ₁^((2σ + 1) * v - 3)

    return add_error(zero(γ), bound)
end

function I_P_dξ_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, normv::Arb, λ::AbstractGLParams{Arb})
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < -1 # Required for integral to converge

    bound =
        C_I_P_dξ(κ, ξ₁, v, λ) *
        normv^(2σ + 1) *
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
    normv::Arb,
    normv_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    return zero(γ) # This is exactly zero
end

function I_P_dγ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    normv::Arb,
    normv_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0 # Required for integral to converge

    bound =
        (2σ + 1) *
        C_I_P(κ, ξ₁, v, λ) *
        normv^2σ *
        normv_dγ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    return add_error(zero(γ), bound)
end

function I_E_dξ_dγ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    normv::Arb,
    normv_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; σ) = λ

    bound = C_I_E_dξ(κ, ξ₁, v, λ) * normv^2σ * normv_dγ * ξ₁^((2σ + 1) * v - 3)

    return add_error(zero(γ), bound)
end

function I_P_dξ_dγ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    normv::Arb,
    normv_dγ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < -1 # Required for integral to converge

    bound =
        (2σ + 1) *
        C_I_P_dξ(κ, ξ₁, v, λ) *
        normv^2σ *
        normv_dγ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    return add_error(zero(γ), bound)
end

function I_E_dκ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    normv::Arb,
    normv_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    return zero(γ) # This is exactly zero
end

function I_P_dκ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    normv::Arb,
    normv_dξ::Arb,
    normv_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ) = λ

    _, _, c = _abc(κ, λ)

    bound =
        (
            C_I_P_dκ_1(κ, ξ₁, v, λ) * normv * ξ₁^(-1) +
            C_I_P_dκ_2(κ, ξ₁, v, λ) * normv +
            C_I_P_dκ_3(κ, ξ₁, v, λ) * normv_dξ +
            C_I_P_dκ_4(κ, ξ₁, v, λ) * normv_dκ * ξ₁^(-1)
        ) *
        normv^2σ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 1)

    return add_error(zero(γ), bound)
end

function I_E_dξ_dκ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    normv::Arb,
    normv_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; σ) = λ

    bound =
        (
            C_I_E_dξ_dκ_1(κ, ξ₁, v, λ) * normv * log(ξ₁) +
            C_I_E_dξ_dκ_2(κ, ξ₁, v, λ) * normv_dκ
        ) *
        normv^2σ *
        ξ₁^((2σ + 1) * v - 3)

    return add_error(zero(γ), bound)
end

function I_P_dξ_dκ_0(
    γ::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    normv::Arb,
    normv_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    (; d, σ, ϵ) = λ

    _, _, c = _abc(κ, λ)

    bound =
        (
            C_I_P_dξ_dκ_1(κ, ξ₁, v, λ) * normv +
            C_I_P_dξ_dκ_2(κ, ξ₁, v, λ) * normv_dκ * ξ₁^(-2)
        ) *
        normv^2σ *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 1)

    return add_error(zero(γ), bound)
end
