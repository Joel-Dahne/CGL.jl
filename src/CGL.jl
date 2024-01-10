module CGL

using Arblib
using ArbExtras
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
include("branch.jl")

end # module CGL
