module CGL

using Arblib
using ArbExtras
using CSV
using DataFrames
using DifferentialEquations
using LinearAlgebra
using NLsolve
using OhMyThreads: tmap, tforeach
using SpecialFunctions
using StaticArrays

import Dates
import Distributed
import ForwardDiff
import IntervalArithmetic
import IntervalArithmetic: BareInterval, bareinterval, inf, sup
import ProgressLogging: @progress, @withprogress, @logprogress

include("CGLBranch/CGLBranch.jl")
include("helper.jl")

include("arb.jl")
include("verify_and_refine_root.jl")
include("special-functions.jl")

include("U.jl")
include("U_expansion.jl")

include("CGLParams.jl")

include("solution_zero/solution.jl")
include("solution_zero/equation.jl")
include("solution_zero/solution_float.jl")
include("solution_zero/solution_taylor.jl")
include("solution_zero/solution_capd.jl")

include("solution_infinity/parameters.jl")
include("solution_infinity/functions.jl")
include("solution_infinity/function_bounds.jl")
include("solution_infinity/norm_bounds_constants.jl")
include("solution_infinity/norm_bounds.jl")
include("solution_infinity/I_bounds.jl")
include("solution_infinity/I.jl")
include("solution_infinity/check_existence.jl")
include("solution_infinity/solution.jl")

include("refine_approximation.jl")
include("refine_approximation_fix_kappa.jl")
include("G.jl")
include("G_solve.jl")

include("branch/branch_points.jl")
include("branch/branch_existence.jl")
include("branch/branch_segment_existence.jl")
include("branch/branch_continuation.jl")
include("branch/data_handling.jl")
include("branch/proof_witness.jl")

include("orchestration/run_branch_points.jl")
include("orchestration/run_branch_points_verification.jl")

using PrecompileTools

@compile_workload begin
    for (j, d) in [(1, 1), (1, 3)]
        μ, γ, κ, ϵ, ξ₁, λ = sverak_params(Float64, j, d)

        G(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)
        G_jacobian_kappa(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)
        G_jacobian_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)

        CGLBranch.branch_epsilon(CGLBranch.sverak_initial(j, d)..., max_steps = 5)

        setprecision(Arb, 128) do
            μ, γ, κ, ϵ, ξ₁, λ = sverak_params(Arb, j, d)

            G(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)
            G_jacobian_kappa(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)
            G_jacobian_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)
        end
    end
end

end # module CGL
