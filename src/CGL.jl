module CGL

using Arblib
using ArbExtras
using CSV
using DataFrames
using DifferentialEquations
using LinearAlgebra
using NLsolve
using SpecialFunctions
using StaticArrays

import Distributed
import ForwardDiff
import IntervalArithmetic
import IntervalArithmetic: BareInterval, bareinterval, inf, sup
import ProgressLogging: @progress, @withprogress, @logprogress

include("CGLBranch/CGLBranch.jl")

include("arb.jl")
include("verify_and_refine_root.jl")
include("special-functions.jl")

include("CGLParams.jl")

include("solution_zero/solution.jl")
include("solution_zero/equation.jl")
include("solution_zero/solution_float.jl")
include("solution_zero/solution_taylor.jl")
include("solution_zero/solution_capd.jl")

include("solution_infinity/functions.jl")
include("solution_infinity/function_bounds.jl")
include("solution_infinity/function_expansions.jl")
include("solution_infinity/solution.jl")
include("solution_infinity/I.jl")
include("solution_infinity/I_bounds.jl")
include("solution_infinity/norm_bounds.jl")
include("solution_infinity/norm_bounds_constants.jl")
include("solution_infinity/check_existence.jl")

include("refine_approximation.jl")
include("G.jl")
include("G_solve.jl")

include("branch/branch_points.jl")
include("branch/branch_existence.jl")
include("branch/branch_segment_existence.jl")
include("branch/branch_continuation.jl")
include("branch/data_handling.jl")

using PrecompileTools

@compile_workload begin
    setprecision(Arb, 128) do
        j, d = 3, 1
        μ, γ, κ, ξ₁, λ = sverak_params(Arb, j, d)

        solution_infinity(γ, κ, ξ₁, λ)
        solution_infinity_jacobian(γ, κ, ξ₁, λ)

        br =
            let λ = CGLBranch.Params(
                    Float64(λ.ϵ),
                    λ.d,
                    Float64(λ.ω),
                    Float64(λ.σ),
                    Float64(λ.δ),
                    30.0,
                )
                CGLBranch.branch(Float64(μ), Float64(κ), λ, max_steps = 5)
            end
    end
end

end # module CGL
