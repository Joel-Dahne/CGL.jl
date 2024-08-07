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
import IntervalArithmetic:
    BareInterval, Interval, bareinterval, interval, inf, sup, isempty_interval, nai
import ProgressLogging: @progress, @withprogress, @logprogress

include("CGLBranch/CGLBranch.jl")

include("arb.jl")
include("interval.jl")
include("helper.jl")
include("special-functions.jl")

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
include("Q_infinity/norm_bounds_constants.jl")
include("Q_infinity/norm_bounds.jl")
include("Q_infinity/I_bounds.jl")
include("Q_infinity/I.jl")
include("Q_infinity/check_existence.jl")
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
include("branch/data_handling.jl")
include("branch/proof_witness.jl")

include("orchestration/run_branch_points.jl")
include("orchestration/run_branch_points_verification.jl")
include("orchestration/run_branch_existence.jl")
include("orchestration/run_branch_continuation.jl")

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
