using ArbExtras
using Arblib
using DifferentialEquations
using CGL
using LinearAlgebra
using SpecialFunctions
using StaticArrays

using Test
using FiniteDifferences

@testset "CGL" verbose = true begin
    include("arb.jl")
    include("verify_and_refine_root.jl")
    include("special-functions.jl")

    include("U.jl")
    include("U_expansion.jl")

    include("solution_zero/solution.jl")
    include("solution_zero/equation.jl")

    include("solution_infinity/parameters.jl")
    include("solution_infinity/functions.jl")
    include("solution_infinity/function_bounds.jl")
    include("solution_infinity/solution.jl")

    include("refine_approximation.jl")
    include("refine_approximation_fix_kappa.jl")
    include("G.jl")
    include("G_solve.jl")
end
