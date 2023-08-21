@testset "solution_infinity" begin
    params = [
        (1.79202 - 1.10431im, 0.493223, gl_params(1, 1.0, 2.3, 0.0, 0.0)),
        (1.71299 - 1.49238im, 0.917383, gl_params(3, 1.0, 1.0, 0.0, 0.0)),
    ]


    ξ₁ = 10.0

    # TODO: Add better tests for this once it is better implemented
    @testset "Parameters $i" for (i, (γ, κ, λ)) in enumerate(params)
        res_F64 = GinzburgLandauSelfSimilarSingular.solution_infinity(γ, κ, ξ₁, λ)
        res_Arb = GinzburgLandauSelfSimilarSingular.solution_infinity(
            Acb(γ),
            Arb(κ),
            Arb(ξ₁),
            gl_params(Arb, λ),
        )

        @test res_F64 ≈ res_Arb

        res_F64_2, res_jacobian_F64 =
            GinzburgLandauSelfSimilarSingular.solution_infinity_jacobian(γ, κ, ξ₁, λ)
        res_Arb_2, res_jacobian_Arb =
            GinzburgLandauSelfSimilarSingular.solution_infinity_jacobian(
                Acb(γ),
                Arb(κ),
                Arb(ξ₁),
                gl_params(Arb, λ),
            )

        @test res_F64_2 ≈ res_Arb_2
        @test res_jacobian_F64 ≈ res_jacobian_Arb
    end
end
