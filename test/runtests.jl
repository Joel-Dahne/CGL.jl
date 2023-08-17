using ArbExtras
using Arblib
using DifferentialEquations
using GinzburgLandauSelfSimilarSingular
using StaticArrays
using Test

@testset "GinzburgLandauSelfSimilarSingular" verbose = true begin
    include("special-functions.jl")

    include("ode_solver/ode_series_solver.jl")

    include("solution_zero/solution.jl")
    include("solution_zero/equation.jl")

    include("solution_infinity/constants.jl")

    include("approximate_parameters.jl")
end
