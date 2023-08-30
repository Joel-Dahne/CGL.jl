"""
    I_E_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, normv::Arb, λ::AbstractGLParams{Arb})

Let
```
I_E(ξ) = ∫_ξ₁^ξ (1 + im * δ) / (1 - im * ϵ) * E(η) * inv(W(η)) * abs(u(η))^2σ * u(η) dη
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
I_P(ξ) = ∫_ξ^∞ (1 + im * δ) / (1 - im * ϵ) * P(η) * inv(W(η)) * abs(u(η))^2σ * u(η) dη
```
This computes a complex ball enclosing `I_P(ξ₁)`.
"""
function I_P_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, normv::Arb, λ::AbstractGLParams{Arb})
    (; d, σ, ϵ) = λ

    c = Acb(0, -κ) / 2Acb(1, -ϵ)
    C_I_P = C_J_P(κ, ξ₁, λ) / abs((2σ + 1) * v - 2 / σ + d - 2)

    bound =
        normv^(2σ + 1) * C_I_P * exp(-real(c) * ξ₁^2) * ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    return add_error(zero(γ), bound)
end

function I_E_dξ_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, normv::Arb, λ::AbstractGLParams{Arb})
    (; σ) = λ

    bound = normv^(2σ + 1) * C_J_E(κ, ξ₁, λ) * ξ₁^((2σ + 1) * v - 3)

    return add_error(zero(γ), bound)
end

function I_P_dξ_0(γ::Acb, κ::Arb, ξ₁::Arb, v::Arb, normv::Arb, λ::AbstractGLParams{Arb})
    (; d, σ, ϵ) = λ

    c = Acb(0, -κ) / 2Acb(1, -ϵ)

    bound =
        normv^(2σ + 1) *
        C_J_P(κ, ξ₁, λ) *
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
    (; d, σ, ϵ) = λ

    c = Acb(0, -κ) / 2Acb(1, -ϵ)
    C_I_P = C_J_P(κ, ξ₁, λ) / abs((2σ + 1) * v - 2 / σ + d - 2)

    bound =
        (2σ + 1) *
        normv^2σ *
        normv_dγ *
        C_I_P *
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

    bound = (2σ + 1) * normv^2σ * normv_dγ * C_J_E(κ, ξ₁, λ) * ξ₁^((2σ + 1) * v - 3)

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
    (; d, σ, ϵ) = λ

    c = Acb(0, -κ) / 2Acb(1, -ϵ)

    bound =
        (2σ + 1) *
        normv^2σ *
        normv_dγ *
        C_J_P(κ, ξ₁, λ) *
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
    normv_dκ::Arb,
    λ::AbstractGLParams{Arb},
)
    # FIXME
    bound = Arb(0)

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
    # FIXME
    bound = Arb(0)

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
    # FIXME
    bound = Arb(0)

    return add_error(zero(γ), bound)
end
