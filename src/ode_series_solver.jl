export ODESeriesSolution, ode_series_solver, ode_series_step

struct ODESeriesSolution{uType,tType}
    u::Vector{uType}
    t::Vector{tType}
end

ODESeriesSolution(u, t) = ODESeriesSolution{eltype(u),eltype(t)}(u, t)

Base.@propagate_inbounds Base.getindex(sol::ODESeriesSolution, i::Union{Integer,Colon}) =
    sol.u[i]

function Base.show(io::IO, sol::ODESeriesSolution{uType,tType}) where {uType,tType}
    println(io, "ODESeriessolution{$uType,$tType}")
    print(io, "$(length(sol.u)) time points")
end

function ode_series_solver(f, u0, tspan, p; degree::Integer = 5, Δt = Arb(0.1))
    Base.require_one_based_indexing(u0)

    t0, tmax = tspan

    u = [u0]
    t = [t0]

    iter = 1
    maxiter = 10000
    while t[end] < tmax
        u0, t0, _ = ode_series_step(f, u0, t0, tmax, p; degree, Δt)

        push!(u, u0)
        push!(t, t0)

        iter += 1
        iter == maxiter && error("too many iterations")
    end

    sol = ODESeriesSolution(u, t)

    return sol
end

function ode_series_step(f, u0, t0, tmax, p; degree::Integer = 5, Δt = Arb(0.1))
    Base.require_one_based_indexing(u0)

    # Compute series at t0
    # TODO: This is written to handle 2nd order ODEs, it might be
    # better to rewrite to handle 1st order ones.
    u_series = f(u0, t0, p; degree)

    # Fix step size
    t = t0 + Δt
    if !(t < tmax)
        t = tmax
        Δt = t - t0
    end

    # Compute value at new point
    # TODO: Don't hard code dimensions
    u = similar(u0)
    u[1], u[3] = Arblib.evaluate2(u_series[1], Δt)
    u[2], u[4] = Arblib.evaluate2(u_series[2], Δt)

    return u, t, u_series
end
