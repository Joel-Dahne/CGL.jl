"""
    _Q_zero_capd(
        Q_ξ₀::SVector{4,BareInterval{Float64}},
        κ::BareInterval{Float64},
        ξ₀::BareInterval{Float64},
        ξ₁::BareInterval{Float64},
        λ::CGLParams{BareInterval{Float64}};
        output_jacobian::Union{Val{false},Val{true}} = Val{false}(),
        jacobian_epsilon::Bool = false,
        tol::Float64 = 1e-11,
    )

Internal function for calling the CAPD program.

It computes the solution to the ODE on the interval ``[0, ξ₁]`` and
returns four real values, the first two are the real and imaginary
values at `ξ₁` and the second two are their derivatives.

If `output_jacobian = Val{true}()` it also computes the Jacobian
w.r.t. `Q_ξ₀` and `κ`. Unless `jacobian_epsilon` is also true, then
it outputs it w.r.t. `Q_ξ₀` and `ϵ`.

- **IMPROVE:** For thin values of `κ` and `ϵ` the performance could be
  improved when computing `Q`. For thin values of `ϵ` it could be
  improved also when computing the Jacobian.
"""
function _Q_zero_capd(
    Q_ξ₀::SVector{4,BareInterval{Float64}},
    κ::BareInterval{Float64},
    ϵ::BareInterval{Float64},
    ξ₀::BareInterval{Float64},
    ξ₁::BareInterval{Float64},
    λ::CGLParams{BareInterval{Float64}};
    output_jacobian::Union{Val{false},Val{true}} = Val{false}(),
    jacobian_epsilon::Bool = false,
    tol::Float64 = 1e-11,
)
    # Build the input to the CAPD program
    input_Q_ξ₀ = ""
    for x in Q_ξ₀
        input_Q_ξ₀ *= "[$(inf(x)), $(sup(x))]\n"
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
    input_tol = "$tol\n"

    input = join([
        input_Q_ξ₀,
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
        Q = fill(IntervalArithmetic.emptyinterval(BareInterval{Float64}), n)
    else
        Q = parse.(
            BareInterval{Float64},
            split(output, "\n"),
        )::Vector{BareInterval{Float64}}
    end

    if output_jacobian isa Val{false}
        return SVector{4,BareInterval{Float64}}(Q)
    else
        return SMatrix{4,5,BareInterval{Float64}}(Q)
    end
end

"""
    Q_zero_capd(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    Q_zero_capd(μ, κ, ϵ, ξ₀, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with four real values, the first two are the real and imaginary
values at `ξ₁` and the second two are their derivatives.

The solution is computed using the rigorous CAPD integrator.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` and CAPD on `[ξ₀, ξ₁]`. For `λ.d != 1` this is
automatically used (`ξ₀ = 1e-2` by default) to handle the removable
singularity at zero.

If the given `ξ₀` gives a non-finite enclosure on `[0, ξ₀]`, then it
tries with half that value. If it fails again it tries to halve it
once more, iterating like this for a maximum of a few times.
"""
Q_zero_capd(μ::Arb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}; tol::Float64 = 1e-11) =
    Q_zero_capd(μ, κ, ϵ, ifelse(isone(λ.d), zero(Arb), Arb(1e-2)), ξ₁, λ; tol)

function Q_zero_capd(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    S = BareInterval{Float64}

    Q_ξ₀ = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        Q_ξ₀ = Q_zero_taylor(μ, κ, ϵ, ξ₀, λ)
        if !all(isfinite, Q_ξ₀)
            iterations = 0
            while !all(isfinite, Q_ξ₀) && iterations < 5
                iterations += 1
                ξ₀ /= 2
                Q_ξ₀ = Q_zero_taylor(μ, κ, ϵ, ξ₀, λ)
            end
        end
        convert(SVector{4,S}, Q_ξ₀)
    else
        SVector{4,S}(μ, bareinterval(0.0), bareinterval(0.0), bareinterval(0.0))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    Q = _Q_zero_capd(
        Q_ξ₀,
        convert(S, κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ);
        tol,
    )

    return Arb.(Q)
end

"""
    Q_zero_jacobian_kappa_capd(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    Q_zero_jacobian_kappa_capd(μ, κ, ϵ, ξ₀, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_zero_capd`](@ref) w.r.t.
the parameters `μ` and `κ`.

Similar to [`Q_zero_capd`](@ref) the solution is computed using the
rigorous CAPD integrator.
"""
Q_zero_jacobian_kappa_capd(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
) = Q_zero_jacobian_kappa_capd(
    μ,
    κ,
    ϵ,
    ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    ξ₁,
    λ;
    tol,
)

function Q_zero_jacobian_kappa_capd(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    S = BareInterval{Float64}

    Q_ξ₀, J_ξ₀ = let
        if !iszero(ξ₀)
            @assert 0 < ξ₀ < ξ₁
            # Integrate system on [0, ξ₀] using Taylor expansion at zero
            Q_ξ₀, J_ξ₀ = Q_zero_jacobian_kappa_taylor(μ, κ, ϵ, ξ₀, λ)
            if !(all(isfinite, Q_ξ₀) && all(isfinite, J_ξ₀))
                iterations = 0
                while !(all(isfinite, Q_ξ₀) && all(isfinite, J_ξ₀)) && iterations < 5
                    iterations += 1
                    ξ₀ /= 2
                    Q_ξ₀, J_ξ₀ = Q_zero_jacobian_kappa_taylor(μ, κ, ϵ, ξ₀, λ)
                end
            end
            Q_ξ₀ = convert(SVector{4,S}, Q_ξ₀)
            J_ξ₀ = convert(SMatrix{4,2,S}, J_ξ₀)
        else
            Q_ξ₀ = SVector{4,S}(μ, bareinterval(0.0), bareinterval(0.0), bareinterval(0.0))
            # Empty integration so the only non-zero derivative is the
            # one of Q_ξ₀[1] w.r.t. μ, which is 1.
            J_ξ₀ = SMatrix{4,2,S}(
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
        # J_ξ₀ now contains derivatives of Q_ξ₀. We want to add a row
        # [0, 1] for the derivative of κ.
        Q_ξ₀, vcat(J_ξ₀, SMatrix{1,2,S}(bareinterval(0.0), bareinterval(1.0)))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    J_ξ₀_ξ₁ = _Q_zero_capd(
        Q_ξ₀,
        convert(S, κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ),
        output_jacobian = Val{true}();
        tol,
    )

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    J = J_ξ₀_ξ₁ * J_ξ₀

    return Arb.(J)
end

"""
    Q_zero_jacobian_epsilon_capd(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    Q_zero_jacobian_epsilon_capd(μ, κ, ϵ, ξ₀, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_zero_capd`](@ref) w.r.t.
the parameters `μ` and `ϵ`.

Similar to [`Q_zero_capd`](@ref) the solution is computed using the
rigorous CAPD integrator.
"""
Q_zero_jacobian_epsilon_capd(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
) = Q_zero_jacobian_epsilon_capd(
    μ,
    κ,
    ϵ,
    ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    ξ₁,
    λ;
    tol,
)

function Q_zero_jacobian_epsilon_capd(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    S = BareInterval{Float64}

    Q_ξ₀, J_ξ₀ = let
        if !iszero(ξ₀)
            @assert 0 < ξ₀ < ξ₁
            # Integrate system on [0, ξ₀] using Taylor expansion at zero
            Q_ξ₀, J_ξ₀ = Q_zero_jacobian_epsilon_taylor(μ, κ, ϵ, ξ₀, λ)
            if !(all(isfinite, Q_ξ₀) && all(isfinite, J_ξ₀))
                iterations = 0
                while !(all(isfinite, Q_ξ₀) && all(isfinite, J_ξ₀)) && iterations < 5
                    iterations += 1
                    ξ₀ /= 2
                    Q_ξ₀, J_ξ₀ = Q_zero_jacobian_kappa_taylor(μ, κ, ϵ, ξ₀, λ)
                end
            end
            Q_ξ₀ = convert(SVector{4,S}, Q_ξ₀)
            J_ξ₀ = convert(SMatrix{4,2,S}, J_ξ₀)
        else
            Q_ξ₀ = SVector{4,S}(μ, bareinterval(0.0), bareinterval(0.0), bareinterval(0.0))
            # Empty integration so the only non-zero derivative is the
            # one of Q_ξ₀[1] w.r.t. μ, which is 1.
            J_ξ₀ = SMatrix{4,2,S}(
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
        # J_ξ₀ now contains derivatives of Q_ξ₀. We want to add a row
        # [0, 1] for the derivative of ϵ.
        Q_ξ₀, vcat(J_ξ₀, SMatrix{1,2,S}(bareinterval(0.0), bareinterval(1.0)))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    J_ξ₀_ξ₁ = _Q_zero_capd(
        Q_ξ₀,
        convert(S, κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ),
        output_jacobian = Val{true}(),
        jacobian_epsilon = true;
        tol,
    )

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    J = J_ξ₀_ξ₁ * J_ξ₀

    return Arb.(J)
end


"""
    _Q_zero_capd_curve(
        Q_ξ₀::SVector{4,BareInterval{Float64}},
        κ::BareInterval{Float64},
        ξ₀::BareInterval{Float64},
        ξ₁::BareInterval{Float64},
        λ::CGLParams{BareInterval{Float64}};
        tol::Float64 = 1e-11,
    )

Internal function for calling the CAPD program.

It computes the solution to the ODE on the interval ``[0, ξ₁]`` and
returns an enclosure of the curve between `ξ₀` and `ξ₁`. It
simultaneously computes the second derivatives.

It returns `ξs, Qs, d2Qs` where
- `ξs::Vector{BareInterval}` contains the `ξ` values
- `Qs::Vector{SVector{4,BareInterval}}` contains the enclosures of the real and
  imaginary parts of `Q` and its derivative for the corresponding `ξ`.
- `d2Qs::Vector{SVector{2,BareInterval}}` contains the enclosurse of the real
  and imaginary parts of the second derivative for the corresponding `ξ`.
"""
function _Q_zero_capd_curve(
    Q_ξ₀::SVector{4,BareInterval{Float64}},
    κ::BareInterval{Float64},
    ϵ::BareInterval{Float64},
    ξ₀::BareInterval{Float64},
    ξ₁::BareInterval{Float64},
    λ::CGLParams{BareInterval{Float64}};
    tol::Float64 = 1e-11,
)
    input_Q_ξ₀ = ""
    for x in Q_ξ₀
        input_Q_ξ₀ *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_params = "$(λ.d)\n"
    for x in [κ, ϵ, λ.ω, λ.σ, λ.δ]
        input_params *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_ξspan = ""
    for x in [ξ₀, ξ₁]
        input_ξspan *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_tol = "$tol\n"

    input = join([input_Q_ξ₀, input_params, input_ξspan, input_tol])

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
        # Return singleton vector with indeterminate enclosure
        ξs = [bareinterval(ξ₀, ξ₁)]
        indet = bareinterval(-Inf, Inf)
        Qs = [SVector(indet, indet, indet, indet)]
        d2Qs = [SVector(indet, indet)]
        abs2_Q_derivative = [indet]
        abs2_Q_derivative2 = [indet]
    else
        res = map(split(output, "\n")) do subinterval
            parse.(BareInterval{Float64}, split(subinterval, ";"))
        end

        ξs = getindex.(res, 1)::Vector{BareInterval{Float64}}
        Qs = [
            SVector(r[2], r[3], r[4], r[5]) for r in res
        ]::Vector{SVector{4,BareInterval{Float64}}}
        d2Qs = [SVector(r[6], r[7]) for r in res]::Vector{SVector{2,BareInterval{Float64}}}
        abs2_Q_derivative = getindex.(res, 8)
        abs2_Q_derivative2 = getindex.(res, 9)
    end

    return ξs, Qs, d2Qs, abs2_Q_derivative, abs2_Q_derivative2
end


"""
    Q_zero_capd_curve(μ, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Similar to [`Q_zero_capd`](@ref) but returns an enclosure of the
solution curve on the entire range, instead of just the value at the
final point. It also returns an enclosure of the second derivatives.

It returns `ξs, Qs, d2Qs` where
- `ξs::Vector{Arb}` contains the `ξ` values.
- `Qs::Vector{SVector{4,Arb}}` contains the enclosures of the real and
  imaginary parts of `Q` and its derivative for the corresponding `ξ`.
- `d2Qs::Vector{SVector{2,Arb}}` contains the enclosurse of the real
  and imaginary parts of the second derivative for the corresponding `ξ`.
"""
Q_zero_capd_curve(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
) = Q_zero_capd_curve(μ, κ, ϵ, ifelse(isone(λ.d), zero(Arb), Arb(1e-2)), ξ₁, λ; tol)

function Q_zero_capd_curve(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    S = BareInterval{Float64}

    Q_ξ₀, d2Q_ξ₀ = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        convert(
            Tuple{SVector{4,S},SVector{2,S}},
            Q_zero_taylor(μ, κ, ϵ, ξ₀, λ, enclose_curve = Val{true}()),
        )
    else
        # d2Q_ξ₀ is not used in this case, so we set it to nai
        SVector{4,S}(
            convert(BareInterval{Float64}, μ),
            bareinterval(0.0),
            bareinterval(0.0),
            bareinterval(0.0),
        ),
        SVector{2,S}(nai(Float64).bareinterval, nai(Float64).bareinterval)
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    ξs, Qs, d2Qs, abs2_Q_derivative, abs2_Q_derivative2 = _Q_zero_capd_curve(
        Q_ξ₀,
        convert(S, κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ);
        tol,
    )

    if !iszero(ξ₀)
        pushfirst!(ξs, bareinterval(0.0, bareinterval(ξ₀)))
        pushfirst!(Qs, Q_ξ₀)
        pushfirst!(d2Qs, d2Q_ξ₀)
        pushfirst!(
            abs2_Q_derivative,
            bareinterval(2) * (Q_ξ₀[3] * Q_ξ₀[1] + Q_ξ₀[4] * Q_ξ₀[2]),
        )
        pushfirst!(
            abs2_Q_derivative2,
            bareinterval(2) * (
                d2Q_ξ₀[1] * Q_ξ₀[1] +
                Q_ξ₀[3]^bareinterval(2) +
                d2Q_ξ₀[2] * Q_ξ₀[2] +
                Q_ξ₀[4]^bareinterval(2)
            ),
        )
    end

    return Arb.(ξs),
    map(Q -> Arb.(Q), Qs),
    map(d2Q -> Arb.(d2Q), d2Qs),
    Arb.(abs2_Q_derivative),
    Arb.(abs2_Q_derivative2)
end
