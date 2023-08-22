@testset "solution_zero" begin
    params = [
        (
            setball(Arb, 0.783077, 1e-9),
            setball(Arb, 0.493223, 1e-9),
            gl_params(Arb, 1, 1.0, 2.3, 0.0, 0.0),
        ),
        (
            setball(Arb, 1.88576, 1e-9),
            setball(Arb, 0.917383, 1e-9),
            gl_params(Arb, 3, 1.0, 1.0, 0.0, 0.0),
        ),
    ]

    ξ₁ = Arb(10.0)

    @testset "Parameters $i" for (i, (μ, κ, λ)) in enumerate(params)
        res_capd = GinzburgLandauSelfSimilarSingular.solution_zero_capd(μ, κ, ξ₁, λ)
        res_float = GinzburgLandauSelfSimilarSingular.solution_zero_float(μ, κ, ξ₁, λ)

        res_capd_2, res_capd_jac =
            GinzburgLandauSelfSimilarSingular.solution_zero_jacobian_capd(μ, κ, ξ₁, λ)
        res_float_2, res_float_jac =
            GinzburgLandauSelfSimilarSingular.solution_zero_jacobian_float(μ, κ, ξ₁, λ)

        # Check that all versions overlap
        @test all(Arblib.overlaps.(res_capd, res_float))
        @test all(Arblib.overlaps.(res_capd, res_capd_2))
        @test all(Arblib.overlaps.(res_capd, res_float_2))
        @test all(Arblib.overlaps.(res_capd_2, res_float))
        @test all(Arblib.overlaps.(res_capd_2, res_float_2))

        @test all(Arblib.overlaps.(res_capd_jac, res_float_jac)) broken = (λ.d != 1)
        if λ.d != 1
            @test res_float_jac \ res_capd_jac ≈ LinearAlgebra.I atol = 1e-2
        end
    end
end
