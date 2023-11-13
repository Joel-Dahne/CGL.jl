module CGL

using Arblib
using ArbExtras
using DifferentialEquations
using LinearAlgebra
using NLsolve
using SpecialFunctions
using StaticArrays

import ForwardDiff
import IntervalArithmetic
import IntervalArithmetic: Interval, interval, inf, sup

include("arb.jl")
include("verify_and_refine_root.jl")
include("special-functions.jl")

include("ode_solver/TaylorModel.jl")
include("ode_solver/ode_series_solver_types.jl")
include("ode_solver/ode_series_solver.jl")

include("GLParams.jl")

include("solution_zero/solution_float.jl")
include("solution_zero/solution_capd.jl")
include("solution_zero/solution.jl")
include("solution_zero/equation.jl")
include("solution_zero/helper.jl")

include("solution_infinity/solution.jl")
include("solution_infinity/I.jl")
include("solution_infinity/I_bounds.jl")
include("solution_infinity/I_second_order.jl")
include("solution_infinity/functions.jl")
include("solution_infinity/function_bounds.jl")
include("solution_infinity/function_expansions.jl")
include("solution_infinity/norm_bounds.jl")
include("solution_infinity/operator_bounds.jl")
include("solution_infinity/check_existence.jl")

include("approximate_parameters.jl")
include("G.jl")

end # module CGL
