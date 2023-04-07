export ODESeriesProblem,
    ODESeriesAutonomusProblem, ODESeriesSecondOrderProblem, ODESeriesSolution

# TODO: Write documentation for the different problem types

abstract type AbstractODESeriesProblem{uType,tType,pType} end

struct ODESeriesProblem{uType,tType,pType} <: AbstractODESeriesProblem{uType,tType,pType}
    f::Any
    u0::Vector{uType}
    tspan::NTuple{2,tType}
    p::pType
end

struct ODESeriesAutonomusProblem{uType,tType,pType} <:
       AbstractODESeriesProblem{uType,tType,pType}
    f::Any
    u0::Vector{uType}
    tspan::NTuple{2,tType}
    p::pType
end

struct ODESeriesSecondOrderProblem{uType,tType,pType} <:
       AbstractODESeriesProblem{uType,tType,pType}
    f::Any
    u0::Vector{NTuple{2,uType}}
    tspan::NTuple{2,tType}
    p::pType
end

ODESeriesProblem(f, u0, tspan, p) =
    ODESeriesProblem{eltype(u0),eltype(tspan),typeof(p)}(f, u0, tspan, p)
ODESeriesAutonomusProblem(f, u0, tspan, p) =
    ODESeriesAutonomusProblem{eltype(u0),eltype(tspan),typeof(p)}(f, u0, tspan, p)
ODESeriesSecondOrderProblem(f, u0, tspan, p) =
    ODESeriesSecondOrderProblem{eltype(u0),eltype(tspan),typeof(p)}(f, u0, tspan, p)

function Base.show(io::IO, prob::ODESeriesProblem{uType,tType}) where {uType,tType}
    println(io, "ODESeriesProblem{$uType,$tType}")
    println(io, "timespan: $(prob.tspan)")
    print(io, "u0: ")
    show(io, prob.u0)
end

function Base.show(io::IO, prob::ODESeriesAutonomusProblem{uType,tType}) where {uType,tType}
    println(io, "ODESeriesAutonomusProblem{$uType,$tType}")
    println(io, "timespan: $(prob.tspan)")
    print(io, "u0: ")
    show(io, prob.u0)
end

function Base.show(
    io::IO,
    prob::ODESeriesSecondOrderProblem{uType,tType},
) where {uType,tType}
    println(io, "ODESeriesSecondOrderProblem{$uType,$tType}")
    println(io, "timespan: $(prob.tspan)")
    print(io, "u0: ")
    show(io, prob.u0)
end

struct ODESeriesSolution{uType,tType}
    u::Vector{uType}
    t::Vector{tType}
end

ODESeriesSolution(u, t) = ODESeriesSolution{eltype(u),eltype(t)}(u, t)

Base.@propagate_inbounds Base.getindex(sol::ODESeriesSolution, i::Union{Integer,Colon}) =
    sol.u[i]

Base.@propagate_inbounds Base.getindex(sol::ODESeriesSolution, ::Colon, i::Integer) =
    getindex.(sol.u, i)

Base.@propagate_inbounds Base.getindex(
    sol::ODESeriesSolution{<:Vector{<:NTuple{2}}},
    ::Colon,
    i::Integer,
    j::Integer,
) = getindex.(getindex.(sol.u, i), j)

function Base.show(io::IO, sol::ODESeriesSolution{uType,tType}) where {uType,tType}
    println(io, "ODESeriesSolution{$uType,$tType}")
    print(io, "$(length(sol.u)) time points")
end
