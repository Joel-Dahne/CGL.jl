"""
    _solve_zero_capd(
        u0::SVector{4,BareInterval{Float64}},
        κ::BareInterval{Float64},
        ξ₀::BareInterval{Float64},
        ξ₁::BareInterval{Float64},
        λ::CGLParams{BareInterval{Float64}},
        output_jacobian::Union{Val{false},Val{true}} = Val{false}(),
    )

Let `u = [a, b, α, β]` be a solution to
[`ivp_zero_real_system`](@ref), but with initial values
```
a(ξ₀) = u0[1]
b(ξ₀) = u0[2]
α(ξ₀) = u0[3]
β(ξ₀) = u0[4]
```
This function computes `u(ξ₁)`.

If `output_jacobian = Val{true}()` it also computes the Jacobian
w.r.t. `[a, b, α, β, κ]`.

- **IMPROVE:** For thin values of `κ` and `λ.ϵ` the performance could
  be improved when computing `u`. For thin values of `λ.ϵ` it could be
  improved also when computing the Jacobian.
"""
function _solve_zero_capd(
    u0::SVector{4,BareInterval{Float64}},
    κ::BareInterval{Float64},
    ξ₀::BareInterval{Float64},
    ξ₁::BareInterval{Float64},
    λ::CGLParams{BareInterval{Float64}},
    output_jacobian::Union{Val{false},Val{true}} = Val{false}(),
)
    input_u0 = ""
    for x in u0
        input_u0 *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_params = "$(λ.d)\n"
    for x in [λ.ω, λ.σ, λ.ϵ, λ.δ, κ]
        input_params *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_ξspan = ""
    for x in [ξ₀, ξ₁]
        input_ξspan *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_output_jacobian = ifelse(output_jacobian isa Val{false}, "0", "1")

    input = join([input_u0, input_params, input_ξspan, input_output_jacobian])

    # IMPROVE: Write directly to stdout of cmd instead of using echo
    program = pkgdir(@__MODULE__, "capd", "build", "ginzburg")
    cmd = pipeline(`echo $input`, `$program`)

    output = readchomp(cmd)

    if contains(output, "Exception")
        n = output_jacobian isa Val{false} ? 4 : 4 + 20
        res = fill(IntervalArithmetic.emptyinterval(BareInterval{Float64}), n)
    else
        res = parse.(
            BareInterval{Float64},
            split(output, "\n"),
        )::Vector{BareInterval{Float64}}
    end

    u = SVector{4,BareInterval{Float64}}(res[1:4])

    if output_jacobian isa Val{false}
        return u
    else
        J = SMatrix{4,5,BareInterval{Float64}}(res[5:end])
        return u, J
    end
end


"""
    solution_zero_capd(μ::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}
    solution_zero_capd(μ::T, κ::T, ξ₀::T, ξ₁::T, λ::CGLParams{T}) where {T}

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes `u(ξ₁)`.

The solution is computed using the rigorous CAPD integrator.

If `p.d != 1` the removable singularity at zero is handled using a
Taylor expansion at zero.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` and CAPD on `[ξ₀, ξ₁]`.
"""
function solution_zero_capd(μ::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}
    if isone(λ.d)
        ξ₀ = zero(ξ₁)
    else
        ξ₀ = convert(T, 1e-2)
    end

    return solution_zero_capd(μ, κ, ξ₀, ξ₁, λ)
end

function solution_zero_capd(μ::T, κ::T, ξ₀::T, ξ₁::T, λ::CGLParams{T}) where {T}
    S = BareInterval{Float64}

    u0 = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        convert(SVector{4,S}, solution_zero_taylor(μ, κ, ξ₀, λ))
    else
        SVector{4,S}(
            convert(BareInterval{Float64}, μ),
            bareinterval(0.0),
            bareinterval(0.0),
            bareinterval(0.0),
        )
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    u = let
        κ = convert(S, κ)
        ξ₀ = convert(S, ξ₀)
        ξ₁ = convert(S, ξ₁)
        λ = CGLParams{S}(λ)

        _solve_zero_capd(u0, κ, ξ₀, ξ₁, λ)
    end

    if T == Float64
        return IntervalArithmetic.mid.(u)
    else
        return convert.(T, u)
    end
end

"""
    solution_zero_jacobian_capd(μ::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}
    solution_zero_jacobian_capd(μ::T, κ::T, ξ₀::T, ξ₁::T, λ::CGLParams{T}) where {T}

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes `u(ξ₁)` as well as the Jacobian w.r.t. `μ` and
`κ`.

The solution is computed using the rigorous CAPD integrator.

If `p.d != 1` the removable singularity at zero is handled using a
Taylor expansion at zero.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` and CAPD on `[ξ₀, ξ₁]`.
"""
function solution_zero_jacobian_capd(μ::T, κ::T, ξ₁::T, λ::CGLParams{T}) where {T}
    if isone(λ.d)
        ξ₀ = zero(ξ₁)
    else
        ξ₀ = convert(T, 1e-2)
    end

    return solution_zero_jacobian_capd(μ, κ, ξ₀, ξ₁, λ)
end

function solution_zero_jacobian_capd(μ::T, κ::T, ξ₀::T, ξ₁::T, λ::CGLParams{T}) where {T}
    S = BareInterval{Float64}

    u0, J1 = let
        if !iszero(ξ₀)
            @assert 0 < ξ₀ < ξ₁
            # Integrate system on [0, ξ₀] using Taylor expansion at zero
            u0, J1 = solution_zero_jacobian_taylor(μ, κ, ξ₀, λ)
            u0 = convert(SVector{4,S}, u0)
            J1 = convert(SMatrix{4,2,S}, J1)
        else
            u0 = SVector{4,S}(
                convert(BareInterval{Float64}, μ),
                bareinterval(0.0),
                bareinterval(0.0),
                bareinterval(0.0),
            )
            # Empty integration so the only non-zero derivative is the
            # one of u0[1] w.r.t. μ, which is 1.
            J1 = SMatrix{4,2,S}(
                bareinterval(1.0),
                bareinterval(0.0),
                bareinterval(0.0),
                bareinterval(0.0),
                bareinterval(0.0),
                bareinterval(0.0),
                bareinterval(0.0),
                bareinterval(0.0),
            )
        end
        # J1 now contains derivatives of [u0[1], u0[2], u0[3], u0[4]].
        # We want to add a row [0, 1] for the derivative of κ.
        u0, vcat(J1, SMatrix{1,2,S}(bareinterval(0.0), bareinterval(1.0)))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    u, J2 = let
        κ = convert(S, κ)
        ξ₀ = convert(S, ξ₀)
        ξ₁ = convert(S, ξ₁)
        λ = CGLParams{S}(λ)

        _solve_zero_capd(u0, κ, ξ₀, ξ₁, λ, Val{true}())
    end

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    J = J2 * J1

    if T == Float64
        return IntervalArithmetic.mid.(u), IntervalArithmetic.mid.(J)
    else
        return convert.(T, u), convert.(T, J)
    end
end
