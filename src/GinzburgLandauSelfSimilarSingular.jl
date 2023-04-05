module GinzburgLandauSelfSimilarSingular

using Arblib
using ArbExtras
using StaticArrays

include("arb.jl")
include("special-functions.jl")

include("ode_series_solver_types.jl")
include("ode_series_solver.jl")

include("GLParams.jl")

include("solution_zero.jl")
include("solution_infinity.jl")

end # module GinzburgLandauSelfSimilarSingular
