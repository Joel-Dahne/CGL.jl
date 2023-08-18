@testset "check_existence" begin
    params = [
        (
            setball(Arb, 0.78308, 1e-10),
            setball(Arb, 0.49323, 1e-10),
            gl_params(Arb, 1, 1.0, 2.3, 0.0, 0.0),
        ),
        (
            setball(Arb, 1.0, 1e-10),
            setball(Arb, 0.45535, 1e-10),
            gl_params(Arb, 3, 1.41727, 1.0, 0.0, 0.0),
        ),
        (
            setball(Arb, 1.0, 1e-10),
            setball(Arb, 0.45535, 1e-10),
            gl_params(Arb, 3, 1.41727, 1.0, 0.01, 0.02),
        ),
    ]

    ξ₁ = Arb(10)
    v = Arb(0.1)

    @testset "Parameters $i" for (i, (μ₀, κ₀, λ)) in enumerate(params)
        if params.d == 1
            sucess1, γ1 = GinzburgLandauSelfSimilarSingular.check_existence_fixed_point(
                μ,
                κ,
                ξ₁,
                v,
                λ,
                non_rigorous = false,
            )

            sucess2, γ2 = GinzburgLandauSelfSimilarSingular.check_existence_fixed_point(
                μ,
                κ,
                ξ₁,
                v,
                λ,
                non_rigorous = true,
            )

            @test sucess1
            @test sucess2
            @test Arblib.overlaps(γ1, γ2)
        else
            @test_throws ArgumentError GinzburgLandauSelfSimilarSingular.check_existence_fixed_point(
                μ,
                κ,
                ξ₁,
                v,
                λ,
            )

            sucess2, γ2 = GinzburgLandauSelfSimilarSingular.check_existence_fixed_point(
                μ,
                κ,
                ξ₁,
                v,
                λ,
                non_rigorous = true,
            )

            @test sucess2
            @test isfinite(γ2)
        end
    end
end
