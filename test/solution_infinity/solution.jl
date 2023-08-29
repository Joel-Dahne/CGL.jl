@testset "solution_infinity" begin
    params = [
        (Acb(1.79202 - 1.10431im), Arb(0.493223), gl_params(Arb, 1, 1.0, 2.3, 0.0, 0.0)),
        (Acb(1.71299 - 1.49238im), Arb(0.917383), gl_params(Arb, 3, 1.0, 1.0, 0.0, 0.0)),
    ]


    ξ₁ = Arb(30.0)

    # TODO: Add better tests for this once it is better implemented

    @testset "solution_infinity" begin
        @testset "Parameters $i" for (i, (γ, κ, λ)) in enumerate(params)
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
        @testset "Parameters $i" for (i, (γ, κ, λ)) in enumerate(params)
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
