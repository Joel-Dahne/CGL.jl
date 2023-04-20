module GinzburgLandauSelfSimilarSingular

using Arblib
using ArbExtras
using SpecialFunctions
using StaticArrays

include("arb.jl")
include("special-functions.jl")
include("TaylorModel.jl")

include("ode_series_solver_types.jl")
include("ode_series_solver.jl")

include("GLParams.jl")

include("solution_zero.jl")
include("solution_zero_helper.jl")

include("solution_infinity.jl")
include("solution_infinity_constants.jl")

end # module GinzburgLandauSelfSimilarSingular
