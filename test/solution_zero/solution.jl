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

        @test all(Arblib.overlaps.(res_capd, res_float))
    end
end
