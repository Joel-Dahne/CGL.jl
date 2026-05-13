using ArbExtras
using Arblib
using CGL
using LinearAlgebra
using SpecialFunctions
using StaticArrays

using Test
using FiniteDifferences

import CGL: P, P_dξ, P_dξ_dξ, P_dξ_dξ_dξ
import CGL: P_dκ, P_dξ_dκ, P_dξ_dξ_dκ
import CGL: P_dϵ, P_dξ_dϵ, P_dξ_dξ_dϵ
import CGL: E, E_dξ, E_dξ_dξ, E_dξ_dξ_dξ
import CGL: E_dκ, E_dξ_dκ
import CGL: E_dϵ, E_dξ_dϵ
import CGL: W
import CGL: J_E, J_E_dξ, J_E_dξ_dξ, J_E_dκ, J_E_dϵ
import CGL: J_P, J_P_dξ, J_P_dξ_dξ, J_P_dκ, J_P_dϵ
import CGL: D, D_dξ, D_dξ_dξ, H, H_dξ, H_dξ_dξ

@testset "CGL" verbose = true begin
    include("arb.jl")
    include("verify_and_refine_root.jl")

    include("U.jl")
    include("U_expansion.jl")

    include("Q_zero/equation.jl")
    include("Q_zero/Q.jl")

    include("Q_infinity/parameters.jl")
    include("Q_infinity/functions.jl")
    include("Q_infinity/function_bounds.jl")
    include("Q_infinity/Q.jl")

    include("refine_approximation.jl")
    include("G.jl")
    include("G_solve.jl")
end
