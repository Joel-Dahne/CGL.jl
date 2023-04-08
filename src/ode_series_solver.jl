export ode_series_solver, ode_series_step

function ode_series_solver(
    prob::AbstractODESeriesProblem;
    degree::Integer = 5,
    Δt::Arb = Arb(0.1),
    skip_remainder::Bool = false,
    verbose::Bool = false,
)
    t0, tmax = prob.tspan

    u = [prob.u0]
    t = [t0]
    if prob isa ODESeriesAutonomusProblem
        β = indeterminate.(prob.u0)
    end

    iter = 1
    maxiter = 10000
    while t[end] < tmax
        if !(t[end] + Δt < tmax)
            Δt = tmax - t[end]
        end

        if prob isa ODESeriesAutonomusProblem
            u_next, t_next, _, β =
                ode_series_step(prob, u[end], t[end], β; degree, Δt, skip_remainder, verbose)
        else
            u_next, t_next, _ = ode_series_step(prob, u[end], t[end]; degree, Δt, skip_remainder, verbose)
        end

        if !isfinite(t_next)
            verbose && @warn "Failed to compute remainder solution at" t[end]
        end

        push!(u, u_next)
        push!(t, t_next)

        iter += 1
        iter == maxiter && error("too many iterations")
    end

    sol = ODESeriesSolution(u, t)

    return sol
end

function ode_series_step(
    prob::ODESeriesSecondOrderProblem,
    u0,
    t0;
    degree::Integer = 5,
    Δt = Arb(0.1),
    skip_remainder = false,
    verbose = false,
)
    # Compute series at t0
    u_expansion = prob.f(u0, t0, prob.p; degree)

    # TODO: Compute remainder
    remainder_coefficient, Δt = let
        [zero(t0) for _ in u_expansion], Δt
    end

    # Construct Taylor model
    I = Arb((t0, t0 + Δt))
    u_model = [TaylorModel(ArbSeries(p, degree = degree + 1), I, t0) for p in u_expansion]
    for i in eachindex(u_model)
        u_model[i].p[end] = remainder_coefficient[i]
    end

    # Compute value at new point
    u = map(p -> Arblib.evaluate2(p.p, Δt), u_model)

    return u, t0 + Δt, u_model
end

function ode_series_step(
    prob::ODESeriesAutonomusProblem,
    u0,
    t0,
    β = indeterminate.(u0);
    degree::Integer = 5,
    Δt = Arb(0.1),
    skip_remainder = false,
    verbose = false,
)
    # Compute expansion at t0
    u_expansion = prob.f(u0, prob.p; degree)

    # Compute remainder, possibly updating Δt
    if skip_remainder
        remainder_coefficient = [zero(t0) for _ in u_expansion]
    else
        remainder_coefficient, Δt = let
            # Enclosure of expansion on [t0, t0 + Δt]
            u_expansion_enclosure = [
                Arb(ArbExtras.extrema_polynomial(p.poly, Arf(0), ubound(Δt))) for
                    p in u_expansion
                    ]

            if !all(isfinite, β)

                # Compute approximate enclosure using an approximate remainder term
                u_tilde_approx = let
                    remainder_coefficient_thin =
                        getindex.(prob.f(u0, prob.p, degree = degree + 1), degree + 1)
                    v =
                        2Δt^(degree + 1) *
                        union.(-remainder_coefficient_thin, remainder_coefficient_thin)

                    u_expansion_enclosure + v
                end

                # Compute guess for upper bound of magnitude for remainder
                # term
                β = abs.(getindex.(prob.f(u_tilde_approx, prob.p, degree = degree + 1), 1))
            end

            for i in eachindex(β)
                if iszero(β[i]) # In practice this probably doesn't happen
                    β[i] = eps(Arb)
                end
            end

            remainder_guess = 2 * Δt^(degree + 1) * union.(-β, β)

            # Compute enclosure given the above guess
            u_tilde = u_expansion_enclosure + remainder_guess

            # Compute remainder given the above enclosure
            remainder_coefficient =
                getindex.(prob.f(u_tilde, prob.p, degree = degree + 1), degree + 1)
            remainder = Δt^(degree + 1) * remainder_coefficient

            # Verify fixed point for guess
            # TODO: Do we need to multiply remainder by Arb((0, 1))?
            if !all(contains.(remainder_guess, remainder))
                # Verification failed, make Δt smaller

                # Want α such that contains.(remainder_guess, α^(degree + 1) * remainder)
                # Heuristically we look at
                # α^(degree + 1) * remainder = remainder_guess
                # and solve for α to get
                # α = (remainder_guess / remainder)^(1 / (degree + 1))
                # In practice we need to take a slightly smaller number.
                # If some part of remainder is exactly zero we also need
                # to skip that.
                min_quotient = minimum(eachindex(remainder)) do i
                    if iszero(abs_ubound(Arb, remainder[i]))
                        Arb(Inf)
                    else
                        abs_ubound(Arb, remainder_guess[i]) / abs_ubound(Arb, remainder[i])
                    end
                end

                α = midpoint(Arb, min_quotient^(1 / (degree + 1)) / 2)

                verbose && @info "Reducing step size" t0 α

                if !all(contains.(remainder_guess, α^(degree + 1) * remainder))
                    q = minimum(abs_ubound.(Arb, remainder_guess) ./ abs_ubound.(Arb, remainder))
                    verbose &&
                        @warn "Verification failed" u_tilde remainder_guess remainder α

                    # If the fixed point is still not okay the return an
                    # indeterminate result
                    remainder_coefficient = indeterminate.(remainder_coefficient)
                end

                Δt *= α
            end

            remainder_coefficient, Δt
        end
    end

    # Construct Taylor model
    I = Arb((t0, t0 + Δt))
    u_model = [TaylorModel(ArbSeries(p, degree = degree + 1), I, t0) for p in u_expansion]
    for i in eachindex(u_model)
        u_model[i].p[end] = remainder_coefficient[i]
    end

    # Compute value at new point
    u = map(p -> p(t0 + Δt), u_model)

    # Compute guess for β for next iteration
    β_next = abs.(remainder_coefficient)

    return u, t0 + Δt, u_model, β_next
end
