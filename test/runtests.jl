using ArbExtras
using Arblib
using DifferentialEquations
using GinzburgLandauSelfSimilarSingular
using StaticArrays
using Test

@testset "GinzburgLandauSelfSimilarSingular" verbose = true begin
    include("special-functions.jl")
    include("solution_zero.jl")
end
