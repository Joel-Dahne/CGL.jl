export ode_series_solver, ode_series_step

function ode_series_solver(
    prob::AbstractODESeriesProblem;
    degree::Integer = 5,
    Δt::Arb = Arb(0.1),
)
    t0, tmax = prob.tspan

    u = [prob.u0]
    t = [t0]

    iter = 1
    maxiter = 10000
    while t[end] < tmax
        u_next, t_next, _ = ode_series_step(prob, u[end], t[end], tmax; degree, Δt)

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
    t0,
    tmax;
    degree::Integer = 5,
    Δt = Arb(0.1),
)
    # Compute series at t0
    u_series = prob.f(u0, t0, prob.p; degree)

    # Fix step size
    t = t0 + Δt
    if !(t < tmax)
        t = tmax
        Δt = t - t0
    end

    # TODO: Compute remainder

    # Compute value at new point. Note the order of the values.
    u = Arblib.evaluate2.(u_series, Δt)

    return u, t, u_series
end

function ode_series_step(
    prob::ODESeriesAutonomusProblem,
    u0,
    t0,
    tmax;
    degree::Integer = 5,
    Δt = Arb(0.1),
)
    # Compute series at t0
    u_series = prob.f(u0, prob.p; degree)

    # Fix step size
    t = t0 + Δt
    if !(t < tmax)
        t = tmax
        Δt = t - t0
    end

    # TODO: Compute remainder

    # Compute value at new point.
    u = [y(Δt) for y in u_series]

    return u, t, u_series
end
