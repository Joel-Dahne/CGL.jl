@testset "enclose_derivative" begin
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

    @testset "Parameters $i" for (i, (μ, κ, λ)) in enumerate(params)
        if λ.d == 1
            res1 = GinzburgLandauSelfSimilarSingular.enclose_derivative_F(
                μ,
                κ,
                ξ₁,
                v,
                λ,
                non_rigorous = false,
            )

            res2 = GinzburgLandauSelfSimilarSingular.enclose_derivative_F(
                μ,
                κ,
                ξ₁,
                v,
                λ,
                non_rigorous = true,
            )

            @test Arblib.overlaps(res1, res2)
        else
            res1 = GinzburgLandauSelfSimilarSingular.enclose_derivative_F(μ, κ, ξ₁, v, λ)

            res2 = GinzburgLandauSelfSimilarSingular.enclose_derivative_F(
                μ,
                κ,
                ξ₁,
                v,
                λ,
                non_rigorous = true,
            )

            @test Arblib.overlaps(res1, res2)
        end
    end
end
