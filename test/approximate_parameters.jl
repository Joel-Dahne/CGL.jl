@testset "approximate_parameters" begin
    params = [
        (0.78308, 0.49323, gl_params(1, 1.0, 2.3, 0.0, 0.0)),
        (1.88590, 0.91742, gl_params(3, 1.0, 1.0, 0.0, 0.0)),
    ]

    ξ₁ = 30.0

    @testset "Parameters $i" for (i, (μ₀, κ₀, λ)) in enumerate(params)
        μ, κ = GinzburgLandauSelfSimilarSingular.approximate_parameters(μ₀, κ₀, ξ₁, λ)

        # Some of the approximations are pretty bad, so we allow for a
        # large error
        @test μ ≈ μ₀ rtol = 1e-4
        @test κ ≈ κ₀ rtol = 1e-4
    end
end
