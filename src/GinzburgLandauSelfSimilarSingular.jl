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
include("special-functions.jl")
include("TaylorModel.jl")

include("ode_series_solver_types.jl")
include("ode_series_solver.jl")

include("GLParams.jl")

include("solution_zero.jl")
include("solution_zero_equation.jl")
include("solution_zero_helper.jl")

include("solution_infinity.jl")
include("solution_infinity_constants.jl")

include("approximate_parameters.jl")
end # module GinzburgLandauSelfSimilarSingular
