module GinzburgLandauSelfSimilarSingular

using Arblib
using ArbExtras
using StaticArrays

include("special-functions.jl")

include("GLParams.jl")

include("solution_zero.jl")
include("solution_infinity.jl")

end # module GinzburgLandauSelfSimilarSingular
