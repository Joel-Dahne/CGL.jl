"""
    _Q_zero_capd(
        Q_ξ₀::SVector{4,Interval{Float64}},
        κ::Interval{Float64},
        ξ₀::Interval{Float64},
        ξ₁::Interval{Float64},
        λ::CGLParams{Interval{Float64}};
        output_jacobian::Union{Val{false},Val{true}} = Val(false),
        jacobian_epsilon::Bool = false,
        output_curve::Union{Val{false},Val{true}} = Val(true),
        tol::Float64 = 1e-11,
    )

Internal function for calling the CAPD program.

It computes the solution to the ODE on the interval ``[ξ₀, ξ₁]`` using
the CAPD C++ library. See `capd/README.md` for more information about
the C++ implementation. The initial value at `ξ₀` is given by `Q_ξ₀`.

By default it returns 4 real values, the first two are the real and
imaginary values at `ξ₁` and the second two are their derivatives.

If `output_jacobian = Val(true)`, then it returns the Jacobian w.r.t.
the initial value `Q_ξ₀` as well as either `κ` or `ϵ` depending on the
argument `wrt_epsilon`.

If `output_curve = Val(true)`, then it returns an enclosure of the
entire curve from `ξ₀` to `ξ₁` as well as values related to the second
derivative. More precisely it returns 5 vectors of the same length:
- `ξs::Vector{Interval}`: Contains intervals in `ξ` covering the
  interval ``[ξ₀, ξ₁]``.
- `Qs::Vector{SVector{4,Interval}}`: Contains enclosures of the real
  and imaginary parts of `Q` and its derivative for the corresponding
  `ξ`.
- `d2Qs::Vector{SVector{2,Interval}}`: Contains enclosures of the real
  and imaginary parts of the second derivative for the corresponding
  `ξ`.
- `abs2_Q_derivative::Vector{Interval}`: Contains an enclosure of the
   derivative of `abs(Q)^2` for the corresponding `ξ`.
- `abs2_Q_derivative2::Vector{Interval}`: Contains an enclosure of the
  second derivative of `abs(Q)^2` for the corresponding `ξ`.
The reason to return the derivatives of `abs(Q)^2` separately is that
these are important for [`count_critical_points`](@ref) and it is
easier to compute accurate enclosures of them directly in the C++
code.
"""
function _Q_zero_capd(
    Q_ξ₀::SVector{4,Interval{Float64}},
    κ::Interval{Float64},
    ϵ::Interval{Float64},
    ξ₀::Interval{Float64},
    ξ₁::Interval{Float64},
    λ::CGLParams{Interval{Float64}};
    output_jacobian::Union{Val{false},Val{true}} = Val(false),
    wrt_epsilon::Bool = false,
    output_curve::Union{Val{false},Val{true}} = Val(false),
    tol::Float64 = 1e-11,
)
    # Only one of these can be set at a time
    @assert !(output_jacobian isa Val{true} && output_curve isa Val{true})

    # Run the C++ program
    exit_success, output = try
        open(`$(pkgdir(@__MODULE__, "capd", "build", "Q_zero"))`, "w+") do io
            # Write initial value
            for x in Q_ξ₀
                println(io, "[$(inf(x)), $(sup(x))]")
            end
            # Write parameters
            println(io, λ.d)
            for x in [κ, ϵ, λ.ω, λ.σ, λ.δ]
                println(io, "[$(inf(x)), $(sup(x))]")
            end
            # Write integration interval
            println(io, "[$(inf(ξ₀)), $(sup(ξ₀))]")
            println(io, "[$(inf(ξ₁)), $(sup(ξ₁))]")
            # Write settings
            println(io, Cint(output_jacobian isa Val{true}))
            println(io, Cint(wrt_epsilon))
            println(io, Cint(output_curve isa Val{true}))
            println(io, tol)
            close(io.in)

            output = readchomp(io)
            exit_success = success(io)
            exit_success, output
        end
    catch e
        e isa ProcessFailedException || rethrow(e)

        # If NaN occurs during the computation the program aborts. We
        # catch this return a failed exit code and empty output.
        false, ""
    end

    if output_jacobian isa Val{true}
        if !exit_success
            J = fill(nai(Float64), 20)
        else
            J = parse.(Interval{Float64}, split(output, "\n"))::Vector{Interval{Float64}}
        end

        return SMatrix{4,5,Interval{Float64}}(J)
    elseif output_curve isa Val{true}
        if !exit_success
            # Return singleton vector with indeterminate enclosure
            ξs = [interval(ξ₀, ξ₁)]
            indet = interval(-Inf, Inf)
            Qs = [SVector(indet, indet, indet, indet)]
            d2Qs = [SVector(indet, indet)]
            abs2_Q_derivative = [indet]
            abs2_Q_derivative2 = [indet]
        else
            res = map(split(output, "\n")) do subinterval
                parse.(Interval{Float64}, split(subinterval, ";"))
            end

            ξs = getindex.(res, 1)::Vector{Interval{Float64}}
            Qs = [
                SVector(r[2], r[3], r[4], r[5]) for r in res
            ]::Vector{SVector{4,Interval{Float64}}}
            d2Qs = [SVector(r[6], r[7]) for r in res]::Vector{SVector{2,Interval{Float64}}}
            abs2_Q_derivative = getindex.(res, 8)
            abs2_Q_derivative2 = getindex.(res, 9)
        end

        return ξs, Qs, d2Qs, abs2_Q_derivative, abs2_Q_derivative2
    else
        if !exit_success
            Q = fill(nai(Float64), 4)
        else
            Q = parse.(Interval{Float64}, split(output, "\n"))::Vector{Interval{Float64}}
        end

        return SVector{4,Interval{Float64}}(Q)
    end
end

"""
    Q_zero_capd(μ, κ, ϵ, ξ₁, λ::CGLParams; ξ₀, tol)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with four real values, the first two are the real and imaginary
values at `ξ₁` and the second two are their derivatives.

The solution is computed using the rigorous CAPD integrator.

If `ξ₀` is non-zero it uses a single Taylor expansion on the interval
`[0, ξ₀]` and the CAPD integrator on `[ξ₀, ξ₁]`. This is needed to
avoid the removable singularity at `ξ = 0` which CAPD cannot handle
directly. For `λ.d = 1` there is no removable singularity and the
default value is `ξ₀ = 0`, otherwise the default value is `ξ₀ = 1e-2`.

If the given `ξ₀` gives a non-finite enclosure on `[0, ξ₀]`, then it
tries with half that value. If it fails again it tries to halve it
once more, iterating like this for a maximum of a few times.
"""
function Q_zero_capd(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    ξ₀::Arb = ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    tol::Float64 = 1e-11,
)
    S = Interval{Float64}

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
            iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
        end
        convert(SVector{4,S}, Q_ξ₀)
    else
        SVector{4,S}(μ, interval(0.0), interval(0.0), interval(0.0))
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
    Q_zero_jacobian_kappa_capd(μ, κ, ϵ, ξ₁, λ::CGLParams; ξ₀, tol)

This function computes the Jacobian of [`Q_zero_capd`](@ref) w.r.t.
the parameters `μ` and `κ`.

In general it works similarly to [`Q_zero_capd`](@ref).
"""
function Q_zero_jacobian_kappa_capd(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    ξ₀::Arb = ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    tol::Float64 = 1e-11,
)
    S = Interval{Float64}

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
                iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
            end
            Q_ξ₀ = convert(SVector{4,S}, Q_ξ₀)
            J_ξ₀ = convert(SMatrix{4,2,S}, J_ξ₀)
        else
            Q_ξ₀ = SVector{4,S}(μ, interval(0.0), interval(0.0), interval(0.0))
            # Empty integration so the only non-zero derivative is the
            # one of Q_ξ₀[1] w.r.t. μ, which is 1.
            J_ξ₀ = SMatrix{4,2,S}(
                interval(1.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
            )
        end
        # J_ξ₀ now contains derivatives of Q_ξ₀. We want to add a row
        # [0, 1] for the derivative of κ.
        Q_ξ₀, vcat(J_ξ₀, SMatrix{1,2,S}(interval(0.0), interval(1.0)))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    J_ξ₀_ξ₁ = _Q_zero_capd(
        Q_ξ₀,
        convert(S, κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ),
        output_jacobian = Val(true);
        tol,
    )

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    J = J_ξ₀_ξ₁ * J_ξ₀

    return Arb.(J)
end

"""
    Q_zero_jacobian_epsilon_capd(μ, κ, ϵ, ξ₁, λ::CGLParams; ξ₀, tol)

This function computes the Jacobian of [`Q_zero_capd`](@ref) w.r.t.
the parameters `μ` and `ϵ`.

In general it works similarly to [`Q_zero_capd`](@ref).
"""
function Q_zero_jacobian_epsilon_capd(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    ξ₀::Arb = ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    tol::Float64 = 1e-11,
)
    S = Interval{Float64}

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
                    Q_ξ₀, J_ξ₀ = Q_zero_jacobian_epsilon_taylor(μ, κ, ϵ, ξ₀, λ)
                end
                iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
            end
            Q_ξ₀ = convert(SVector{4,S}, Q_ξ₀)
            J_ξ₀ = convert(SMatrix{4,2,S}, J_ξ₀)
        else
            Q_ξ₀ = SVector{4,S}(μ, interval(0.0), interval(0.0), interval(0.0))
            # Empty integration so the only non-zero derivative is the
            # one of Q_ξ₀[1] w.r.t. μ, which is 1.
            J_ξ₀ = SMatrix{4,2,S}(
                interval(1.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
            )
        end
        # J_ξ₀ now contains derivatives of Q_ξ₀. We want to add a row
        # [0, 1] for the derivative of ϵ.
        Q_ξ₀, vcat(J_ξ₀, SMatrix{1,2,S}(interval(0.0), interval(1.0)))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    J_ξ₀_ξ₁ = _Q_zero_capd(
        Q_ξ₀,
        convert(S, κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ),
        output_jacobian = Val(true),
        wrt_epsilon = true;
        tol,
    )

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    J = J_ξ₀_ξ₁ * J_ξ₀

    return Arb.(J)
end

"""
    Q_zero_capd_curve(μ, κ, ϵ, ξ₁, λ::CGLParams; ξ₀, tol)

Similar to [`Q_zero_capd`](@ref) but returns an enclosure of the
solution curve on the entire range, instead of just the value at the
final point. It also returns values related to the second derivative
and `abs(Q)^2`, since these are needed in
[`count_critical_points`](@ref).

It returns 5 vectors, all of the same length:
- `ξs::Vector{Interval}`: Contains intervals in `ξ` covering the
  interval ``[ξ₀, ξ₁]``.
- `Qs::Vector{SVector{4,Interval}}`: Contains enclosures of the real
  and imaginary parts of `Q` and its derivative for the corresponding
  `ξ`.
- `d2Qs::Vector{SVector{2,Interval}}`: Contains enclosures of the real
  and imaginary parts of the second derivative for the corresponding
  `ξ`.
- `abs2_Q_derivative::Vector{Interval}`: Contains an enclosure of the
   derivative of `abs(Q)^2` for the corresponding `ξ`.
- `abs2_Q_derivative2::Vector{Interval}`: Contains an enclosure of the
  second derivative of `abs(Q)^2` for the corresponding `ξ`.

In general it works similarly to [`Q_zero_capd`](@ref).
"""
function Q_zero_capd_curve(
    μ::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    ξ₀::Arb = ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    tol::Float64 = 1e-11,
)
    S = Interval{Float64}

    Q_ξ₀, d2Q_ξ₀ = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        convert(
            Tuple{SVector{4,S},SVector{2,S}},
            Q_zero_taylor(μ, κ, ϵ, ξ₀, λ, enclose_curve = Val(true)),
        )
    else
        # d2Q_ξ₀ is not used in this case, so we set it to nai
        SVector{4,S}(convert(S, μ), interval(0.0), interval(0.0), interval(0.0)),
        SVector{2,S}(nai(Float64), nai(Float64))
    end

    # Integrate system on [ξ₀, ξ₁] using capd
    ξs, Qs, d2Qs, abs2_Q_derivative, abs2_Q_derivative2 = _Q_zero_capd(
        Q_ξ₀,
        convert(S, κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ),
        output_curve = Val(true);
        tol,
    )

    if !iszero(ξ₀)
        pushfirst!(ξs, interval(0.0, interval(ξ₀)))
        pushfirst!(Qs, Q_ξ₀)
        pushfirst!(d2Qs, d2Q_ξ₀)
        pushfirst!(abs2_Q_derivative, interval(2) * (Q_ξ₀[3] * Q_ξ₀[1] + Q_ξ₀[4] * Q_ξ₀[2]))
        pushfirst!(
            abs2_Q_derivative2,
            interval(2) * (
                d2Q_ξ₀[1] * Q_ξ₀[1] +
                Q_ξ₀[3]^interval(2) +
                d2Q_ξ₀[2] * Q_ξ₀[2] +
                Q_ξ₀[4]^interval(2)
            ),
        )
    end

    return Arb.(ξs),
    map(Q -> Arb.(Q), Qs),
    map(d2Q -> Arb.(d2Q), d2Qs),
    Arb.(abs2_Q_derivative),
    Arb.(abs2_Q_derivative2)
end
