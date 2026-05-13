module CGL

using Arblib
using ArbExtras
using CSV
using DataFrames
using LinearAlgebra
using NonlinearSolve
using OhMyThreads: tmap, tforeach
using OrdinaryDiffEqVerner
using SpecialFunctions
using StaticArrays

import Dates
import Distributed
import ForwardDiff
import IntervalArithmetic: Interval, interval, nai, inf, sup
import ProgressLogging: @progress, @withprogress, @logprogress

include("CGLBranch/CGLBranch.jl")

include("assert_proof.jl")

include("arb.jl")
include("helper.jl")

include("verify_and_refine_root.jl")

include("U.jl")
include("U_expansion.jl")

include("CGLParams.jl")

include("Q_zero/equation.jl")
include("Q_zero/Q.jl")
include("Q_zero/Q_float.jl")
include("Q_zero/Q_taylor.jl")
include("Q_zero/Q_capd.jl")

include("Q_infinity/parameters.jl")
include("Q_infinity/functions.jl")
include("Q_infinity/function_bounds.jl")
include("Q_infinity/I_bounds.jl")
include("Q_infinity/norm_bounds_constants.jl")
include("Q_infinity/norm_bounds.jl")
include("Q_infinity/I_P.jl")
include("Q_infinity/Q.jl")
include("Q_infinity/verify_monotonicity.jl")

include("refine_approximation.jl")
include("G.jl")
include("G_solve.jl")
include("count_critical_points.jl")

include("branch/branch_points.jl")
include("branch/branch_existence.jl")
include("branch/branch_segment_existence.jl")
include("branch/branch_continuation.jl")
include("branch/branch_critical_points.jl")
include("branch/data_handling.jl")
include("branch/proof_witness.jl")

include("orchestration/run_branch_points.jl")
include("orchestration/run_branch_points_verification.jl")
include("orchestration/run_branch_existence.jl")
include("orchestration/run_branch_continuation.jl")
include("orchestration/run_branch_critical_points.jl")

using PrecompileTools

@compile_workload begin
    for (j, d) in [(1, 1), (1, 3)]
        μ, γ, κ, ϵ, ξ₁, Λ = sverak_params(Float64, j, d)

        G(μ, real(γ), imag(γ), κ, ϵ, ξ₁, Λ)
        G_jacobian_kappa(μ, real(γ), imag(γ), κ, ϵ, ξ₁, Λ)
        G_jacobian_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁, Λ)

        CGLBranch.branch_epsilon(CGLBranch.sverak_initial(j, d)..., max_steps = 5)

        setprecision(Arb, 128) do
            μ, γ, κ, ϵ, ξ₁, Λ = sverak_params(Arb, j, d)

            G(μ, real(γ), imag(γ), κ, ϵ, ξ₁, Λ)
            G_jacobian_kappa(μ, real(γ), imag(γ), κ, ϵ, ξ₁, Λ)
            G_jacobian_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁, Λ)
        end
    end
end

end # module CGL
