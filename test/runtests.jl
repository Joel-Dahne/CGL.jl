using ArbExtras
using Arblib
using DifferentialEquations
using GinzburgLandauSelfSimilarSingular
using StaticArrays
using Test

@testset "GinzburgLandauSelfSimilarSingular" verbose = true begin
    include("special-functions.jl")
    include("ode_series_solver.jl")
    include("solution_zero_equation.jl")
    include("solution_zero.jl")
    include("solution_infinity_constants.jl")
    include("approximate_parameters.jl")
end
