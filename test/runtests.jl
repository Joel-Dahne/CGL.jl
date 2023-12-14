using ArbExtras
using Arblib
using DifferentialEquations
using CGL
using LinearAlgebra
using StaticArrays

using Test
using FiniteDifferences

@testset "CGL" verbose = true begin
    include("arb.jl")
    include("special-functions.jl")
    include("verify_and_refine_root.jl")

    include("ode_solver/ode_series_solver.jl")

    include("solution_zero/solution.jl")
    include("solution_zero/equation.jl")

    include("solution_infinity/solution.jl")
    include("solution_infinity/functions.jl")
    include("solution_infinity/function_expansions.jl")
    include("solution_infinity/function_bounds.jl")
    include("solution_infinity/constants.jl")

    include("approximate_parameters.jl")
    include("G.jl")
    include("G_solve.jl")
end
