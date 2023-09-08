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

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0 # Required for integral to converge

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

    @assert (2σ + 1) * v - 2 / σ + d - 3 < -1 # Required for integral to converge

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

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0 # Required for integral to converge

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

    @assert (2σ + 1) * v - 2 / σ + d - 3 < -1 # Required for integral to converge

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
    (; d, σ, ϵ) = λ

    c = Acb(0, -κ) / 2Acb(1, -ϵ)
    C_I_P = C_J_P(κ, ξ₁, λ) / abs((2σ + 1) * v - 2 / σ + d - 2)

    @assert (2σ + 1) * v - 2 / σ + d - 2 < 0 # Required for integral to converge

    bound_1 =
        (2σ + 1) *
        normv^2σ *
        normv_dκ *
        C_I_P *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 2)

    # FIXME
    #@assert XXX < 0 # Required for integral to converge

    bound_2 = normv^(2σ + 1) * C_J_P_dκ(κ, ξ₁, λ) * 0 # TODO

    bound = bound_1 + bound_2

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

    bound_1 = (2σ + 1) * normv^2σ * normv_dκ * C_J_E(κ, ξ₁, λ) * ξ₁^((2σ + 1) * v - 3)

    bound_2 = normv^(2σ + 1) * C_J_E_dκ(κ, ξ₁, λ) * log(ξ₁) * ξ₁^((2σ + 1) * v - 3)

    bound = bound_1 + bound_2

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

    c = Acb(0, -κ) / 2Acb(1, -ϵ)

    @assert (2σ + 1) * v - 2 / σ + d - 3 < -1 # Required for integral to converge

    bound_1 =
        (2σ + 1) *
        normv^2σ *
        normv_dκ *
        C_J_P(κ, ξ₁, λ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 3)

    # FIXME
    #@assert (2σ + 1) * v - 2 / σ + d - 1 < -1 # Required for integral to converge

    bound_2 =
        normv^(2σ + 1) *
        C_J_E_dκ(κ, ξ₁, λ) *
        exp(-real(c) * ξ₁^2) *
        ξ₁^((2σ + 1) * v - 2 / σ + d - 1)

    bound = bound_1 + bound_2

    return add_error(zero(γ), bound)
end
