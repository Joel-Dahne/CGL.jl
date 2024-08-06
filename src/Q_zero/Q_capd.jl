"""
    _solve_zero_capd_curve(
        u0::SVector{4,BareInterval{Float64}},
        κ::BareInterval{Float64},
        ξ₀::BareInterval{Float64},
        ξ₁::BareInterval{Float64},
        λ::CGLParams{BareInterval{Float64}};
        tol::Float64 = 1e-11,
    )

Let `u = [a, b, α, β]` be a solution to
[`ivp_zero_real_system`](@ref), but with initial values
```
a(ξ₀) = u0[1]
b(ξ₀) = u0[2]
α(ξ₀) = u0[3]
β(ξ₀) = u0[4]
```
This function computes an enclosure of the curve between `ξ₀` and
`ξ₁`. It simultaneously computes the second derivatives of `a` and `b`
(i.e. the first derivatives of `α` and `β`).

It returns `ξs, us, d2us` where
- `ξs::Vector{BareInterval}` contains the `ξ` values
- `us::Vector{SVector{4,BareInterval}}` contains the enclosures of `a,
  b, α, β` for the corresponding `ξ`
- `du2s::Vector{SVector{2,BareInterval}}` contains the enclosure of the
  derivative of `α` and `β` for the corresponding `ξ`.
"""
function _solve_zero_capd_curve(
    u0::SVector{4,BareInterval{Float64}},
    κ::BareInterval{Float64},
    ϵ::BareInterval{Float64},
    ξ₀::BareInterval{Float64},
    ξ₁::BareInterval{Float64},
    λ::CGLParams{BareInterval{Float64}};
    tol::Float64 = 1e-11,
)
    input_u0 = ""
    for x in u0
        input_u0 *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_params = "$(λ.d)\n"
    for x in [κ, ϵ, λ.ω, λ.σ, λ.δ]
        input_params *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_ξspan = ""
    for x in [ξ₀, ξ₁]
        input_ξspan *= "[$(inf(x)), $(sup(x))]\n"
    end
    # IMPROVE: Consider choosing the tolerance depending on the input
    input_tol = "$tol\n"

    input = join([input_u0, input_params, input_ξspan, input_tol])

    # IMPROVE: Write directly to stdout of cmd instead of using echo
    program = pkgdir(@__MODULE__, "capd", "build", "ginzburg-curve")
    cmd = pipeline(`echo $input`, `$program`)

    output = try
        readchomp(cmd)
    catch e
        # If NaN occurs during the computation the program aborts. We
        # catch this and handle it in the same way as if an exception
        # was thrown during the computations.
        e isa ProcessFailedException || rethrow(e)

        "Exception"
    end

    if contains(output, "Exception")
        ξs = BareInterval{Float64}[]
        us = SVector{4,BareInterval{Float64}}[]
    else
        res = map(split(output, "\n")) do subinterval
            parse.(BareInterval{Float64}, split(subinterval, ";"))
        end

        ξs = getindex.(res, 1)::Vector{BareInterval{Float64}}
        us = [
            SVector(r[2], r[3], r[4], r[5]) for r in res
        ]::Vector{SVector{4,BareInterval{Float64}}}
        d2us = [SVector(r[6], r[7]) for r in res]::Vector{SVector{2,BareInterval{Float64}}}

    end

    return ξs, us, d2us
end


"""
    _solve_zero_capd(
        u0::SVector{4,BareInterval{Float64}},
        κ::BareInterval{Float64},
        ξ₀::BareInterval{Float64},
        ξ₁::BareInterval{Float64},
        λ::CGLParams{BareInterval{Float64}};
        output_jacobian::Union{Val{false},Val{true}} = Val{false}(),
        jacobian_epsilon::Bool = false,
        tol::Float64 = 1e-11,
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
w.r.t. `[a, b, α, β, κ]`. Unless `jacobian_epsilon` is also true, then
it outputs it w.r.t. `[a, b, α, β, ϵ]`.

- **IMPROVE:** For thin values of `κ` and `λ.ϵ` the performance could
  be improved when computing `u`. For thin values of `λ.ϵ` it could be
  improved also when computing the Jacobian.
"""
function _solve_zero_capd(
    u0::SVector{4,BareInterval{Float64}},
    κ::BareInterval{Float64},
    ϵ::BareInterval{Float64},
    ξ₀::BareInterval{Float64},
    ξ₁::BareInterval{Float64},
    λ::CGLParams{BareInterval{Float64}};
    output_jacobian::Union{Val{false},Val{true}} = Val{false}(),
    jacobian_epsilon::Bool = false,
    tol::Float64 = 1e-11,
)
    input_u0 = ""
    for x in u0
        input_u0 *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_params = "$(λ.d)\n"
    for x in [κ, ϵ, λ.ω, λ.σ, λ.δ]
        input_params *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_ξspan = ""
    for x in [ξ₀, ξ₁]
        input_ξspan *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_output_jacobian = ifelse(output_jacobian isa Val{true}, "1\n", "0\n")
    input_jacobian_epsilon = ifelse(jacobian_epsilon, "1\n", "0\n")
    # IMPROVE: Consider choosing the tolerance depending on the input
    input_tol = "$tol\n"

    input = join([
        input_u0,
        input_params,
        input_ξspan,
        input_output_jacobian,
        input_jacobian_epsilon,
        input_tol,
    ])

    # IMPROVE: Write directly to stdout of cmd instead of using echo
    program = pkgdir(@__MODULE__, "capd", "build", "ginzburg")
    cmd = pipeline(`echo $input`, `$program`)

    output = try
        readchomp(cmd)
    catch e
        # If NaN occurs during the computation the program aborts. We
        # catch this and handle it in the same way as if an exception
        # was thrown during the computations.
        e isa ProcessFailedException || rethrow(e)

        "Exception"
    end

    if contains(output, "Exception")
        n = output_jacobian isa Val{false} ? 4 : 20
        res = fill(IntervalArithmetic.emptyinterval(BareInterval{Float64}), n)
    else
        res = parse.(
            BareInterval{Float64},
            split(output, "\n"),
        )::Vector{BareInterval{Float64}}
    end

    if output_jacobian isa Val{false}
        return SVector{4,BareInterval{Float64}}(res)
    else
        return SMatrix{4,5,BareInterval{Float64}}(res)
    end
end

"""
    Q_zero_capd_curve(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Similar to [`Q_zero_capd`](@ref) but returns an enclosure of the
solution curve on the entire range, instead of just the value at the
final point. It also reeturns an enclosure of the second derivative of
`a` and `b`.

It returns `ξs, us, d2us` where
- `ξs::Vector{BareInterval}` contains the `ξ` values
- `us::Vector{SVector{4,BareInterval}}` contains the enclosures of `a,
  b, α, β` for the corresponding `ξ`
- `du2s::Vector{SVector{2,BareInterval}}` contains the enclosure of the
  derivative of `α` and `β` for the corresponding `ξ`.
"""
function Q_zero_capd_curve(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    if isone(λ.d)
        ξ₀ = zero(ξ₁)
    else
        ξ₀ = convert(T, 1e-2)
    end

    return Q_zero_capd_curve(μ, κ, ϵ, ξ₀, ξ₁, λ; tol)
end

function Q_zero_capd_curve(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₀::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    S = BareInterval{Float64}

    ξ0 = bareinterval(0.0, convert(BareInterval{Float64}, ξ₀))

    u0, d2u0 = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        convert(
            Tuple{SVector{4,S},SVector{2,S}},
            Q_zero_taylor(μ, κ, ϵ, ξ₀, λ, enclose_curve = Val{true}()),
        )
    else
        # d2u0 is not used in this case, so we set it to an entire
        # interval
        SVector{4,S}(
            convert(BareInterval{Float64}, μ),
            bareinterval(0.0),
            bareinterval(0.0),
            bareinterval(0.0),
        ),
        SVector{2,S}(bareinterval(-Inf, Inf), bareinterval(-Inf, Inf))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    ξs, us, d2us =
        let κ = convert(S, κ),
            ϵ = convert(S, ϵ),
            ξ₀ = convert(S, ξ₀),
            ξ₁ = convert(S, ξ₁),
            λ = CGLParams{S}(λ)

            _solve_zero_capd_curve(u0, κ, ϵ, ξ₀, ξ₁, λ; tol)
        end

    if !iszero(ξ₀)
        pushfirst!(ξs, ξ0)
        pushfirst!(us, u0)
        pushfirst!(d2us, d2u0)
    end

    if T == Float64
        return IntervalArithmetic.mid.(ξs),
        map(u -> IntervalArithmetic.mid.(u), us),
        map(u -> IntervalArithmetic.mid.(u), d2us)
    else
        return convert.(T, ξs), map(u -> convert.(T, u), us), map(u -> convert.(T, u), d2us)
    end
end


"""
    Q_zero_capd(μ::T, κ::T, ϵ::T, ξ₁::T, λ::CGLParams{T}; tol::Float64 = 1e-11) where {T}
    Q_zero_capd(μ::T, κ::T, ϵ::T, ξ₀::T, ξ₁::T, λ::CGLParams{T}; tol::Float64 = 1e-11) where {T}

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes `u(ξ₁)`.

The solution is computed using the rigorous CAPD integrator.

If `p.d != 1` the removable singularity at zero is handled using a
Taylor expansion at zero.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` and CAPD on `[ξ₀, ξ₁]`.
"""
function Q_zero_capd(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    if isone(λ.d)
        ξ₀ = zero(ξ₁)
    else
        ξ₀ = convert(T, 1e-2)
    end

    return Q_zero_capd(μ, κ, ϵ, ξ₀, ξ₁, λ; tol)
end

function Q_zero_capd(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₀::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    S = BareInterval{Float64}

    u0 = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        convert(SVector{4,S}, Q_zero_taylor(μ, κ, ϵ, ξ₀, λ))
    else
        SVector{4,S}(
            convert(BareInterval{Float64}, μ),
            bareinterval(0.0),
            bareinterval(0.0),
            bareinterval(0.0),
        )
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    u =
        let κ = convert(S, κ),
            ϵ = convert(S, ϵ),
            ξ₀ = convert(S, ξ₀),
            ξ₁ = convert(S, ξ₁),
            λ = CGLParams{S}(λ)

            _solve_zero_capd(u0, κ, ϵ, ξ₀, ξ₁, λ; tol)
        end

    if T == Float64
        return IntervalArithmetic.mid.(u)
    else
        return convert.(T, u)
    end
end

"""
    Q_zero_jacobian_kappa_capd(μ::T, κ::T, ϵ::T, ξ₁::T, λ::CGLParams{T}; tol::Float64 = 1e-11) where {T}
    Q_zero_jacobian_kappa_capd(μ::T, κ::T, ϵ::T, ξ₀::T, ξ₁::T, λ::CGLParams{T}; tol::Float64 = 1e-11) where {T}

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes the Jacobian w.r.t. `μ` and `κ`.

The solution is computed using the rigorous CAPD integrator.

If `p.d != 1` the removable singularity at zero is handled using a
Taylor expansion at zero.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` and CAPD on `[ξ₀, ξ₁]`.
"""
function Q_zero_jacobian_kappa_capd(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    if isone(λ.d)
        ξ₀ = zero(ξ₁)
    else
        ξ₀ = convert(T, 1e-2)
    end

    return Q_zero_jacobian_kappa_capd(μ, κ, ϵ, ξ₀, ξ₁, λ; tol)
end

function Q_zero_jacobian_kappa_capd(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₀::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    S = BareInterval{Float64}

    u0, J1 = let
        if !iszero(ξ₀)
            @assert 0 < ξ₀ < ξ₁
            # Integrate system on [0, ξ₀] using Taylor expansion at zero
            u0, J1 = Q_zero_jacobian_kappa_taylor(μ, κ, ϵ, ξ₀, λ)
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
    J2 =
        let κ = convert(S, κ),
            ϵ = convert(S, ϵ),
            ξ₀ = convert(S, ξ₀),
            ξ₁ = convert(S, ξ₁),
            λ = CGLParams{S}(λ)

            _solve_zero_capd(u0, κ, ϵ, ξ₀, ξ₁, λ, output_jacobian = Val{true}(); tol)
        end

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    J = J2 * J1

    if T == Float64
        return IntervalArithmetic.mid.(J)
    else
        return convert.(T, J)
    end
end

"""
    Q_zero_jacobian_epsilon_capd(μ::T, κ::T, ϵ::T, ξ₁::T, λ::CGLParams{T}; tol::Float64 = 1e-11) where {T}
    Q_zero_jacobian_epsilon_capd(μ::T, κ::T, ϵ::T, ξ₀::T, ξ₁::T, λ::CGLParams{T}; tol::Float64 = 1e-11) where {T}

Let `u = [a, b, α, β]` be a solution to [`ivp_zero_real_system`](@ref)
This function computes `u(ξ₁)` as well as the Jacobian w.r.t. `μ` and
`ϵ`.

The solution is computed using the rigorous CAPD integrator.

If `p.d != 1` the removable singularity at zero is handled using a
Taylor expansion at zero.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` and CAPD on `[ξ₀, ξ₁]`.
"""
function Q_zero_jacobian_epsilon_capd(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    if isone(λ.d)
        ξ₀ = zero(ξ₁)
    else
        ξ₀ = convert(T, 1e-2)
    end

    return Q_zero_jacobian_epsilon_capd(μ, κ, ϵ, ξ₀, ξ₁, λ; tol)
end

function Q_zero_jacobian_epsilon_capd(
    μ::T,
    κ::T,
    ϵ::T,
    ξ₀::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    S = BareInterval{Float64}

    u0, J1 = let
        if !iszero(ξ₀)
            @assert 0 < ξ₀ < ξ₁
            # Integrate system on [0, ξ₀] using Taylor expansion at zero
            u0, J1 = Q_zero_jacobian_epsilon_taylor(μ, κ, ϵ, ξ₀, λ)
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
        # We want to add a row [0, 1] for the derivative of ϵ.
        u0, vcat(J1, SMatrix{1,2,S}(bareinterval(0.0), bareinterval(1.0)))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    J2 =
        let κ = convert(S, κ),
            ϵ = convert(S, ϵ),
            ξ₀ = convert(S, ξ₀),
            ξ₁ = convert(S, ξ₁),
            λ = CGLParams{S}(λ)

            _solve_zero_capd(
                u0,
                κ,
                ϵ,
                ξ₀,
                ξ₁,
                λ,
                output_jacobian = Val{true}(),
                jacobian_epsilon = true;
                tol,
            )
        end

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    J = J2 * J1

    if T == Float64
        return IntervalArithmetic.mid.(J)
    else
        return convert.(T, J)
    end
end
