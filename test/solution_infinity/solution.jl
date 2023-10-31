@testset "solution_infinity" begin
    params = [GinzburgLandauSelfSimilarSingular._params.(Arb, 1, d) for d in [1, 3]]

    @testset "solution_infinity" begin
        @testset "Parameters $i" for (i, (μ, γ, κ, ξ₁, λ)) in enumerate(params)
            res_F64 = GinzburgLandauSelfSimilarSingular.solution_infinity(
                Complex{Float64}(γ),
                Float64(κ),
                Float64(ξ₁),
                gl_params(Float64, λ),
            )
            res_Arb = GinzburgLandauSelfSimilarSingular.solution_infinity(γ, κ, ξ₁, λ)

            @test res_F64 ≈ ComplexF64.(res_Arb)
        end
    end

    @testset "solution_infinity_jacobian" begin
        @testset "Parameters $i" for (i, (μ, γ, κ, ξ₁, λ)) in enumerate(params)
            res_jacobian_F64 = GinzburgLandauSelfSimilarSingular.solution_infinity_jacobian(
                Complex{Float64}(γ),
                Float64(κ),
                Float64(ξ₁),
                gl_params(Float64, λ),
            )
            res_jacobian_Arb =
                GinzburgLandauSelfSimilarSingular.solution_infinity_jacobian(γ, κ, ξ₁, λ)

            @test res_jacobian_F64 ≈ ComplexF64.(res_jacobian_Arb)
        end
    end
end
