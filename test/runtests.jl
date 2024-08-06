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

    include("Q_zero/equation.jl")
    include("Q_zero/Q.jl")

    include("Q_infinity/parameters.jl")
    include("Q_infinity/functions.jl")
    include("Q_infinity/function_bounds.jl")
    include("Q_infinity/Q.jl")

    include("refine_approximation.jl")
    include("G.jl")
    include("G_solve.jl")
end
