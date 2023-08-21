module GinzburgLandauSelfSimilarSingular

using Arblib
using ArbExtras
using DifferentialEquations
using NLsolve
using SpecialFunctions
using StaticArrays

import IntervalArithmetic
import IntervalArithmetic: Interval

include("arb.jl")
include("verify_and_refine_root.jl")
include("special-functions.jl")

include("ode_solver/TaylorModel.jl")
include("ode_solver/ode_series_solver_types.jl")
include("ode_solver/ode_series_solver.jl")

include("GLParams.jl")

include("solution_zero/solution.jl")
include("solution_zero/equation.jl")
include("solution_zero/helper.jl")

include("solution_infinity/solution.jl")
include("solution_infinity/functions.jl")
include("solution_infinity/function_bounds.jl")
include("solution_infinity/operator_bounds.jl")
include("solution_infinity/check_existence.jl")

include("approximate_parameters.jl")
include("G.jl")
include("enclose_derivative.jl")

end # module GinzburgLandauSelfSimilarSingular
