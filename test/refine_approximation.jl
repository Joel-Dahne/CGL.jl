@testset "refine_approximation" begin
    params = [
        (0.78308, 0.49323, CGLParams(1, 1.0, 2.3, 0.0, 0.0)),
        (1.88590, 0.91742, CGLParams(3, 1.0, 1.0, 0.0, 0.0)),
    ]

    ξ₁ = 30.0

    @testset "Parameters $i" for (i, (μ₀, κ₀, λ)) in enumerate(params)
        μ, γ, κ = CGL.refine_approximation(μ₀, κ₀, ξ₁, λ)

        @test μ ≈ μ₀ rtol = 1e-3
        @test κ ≈ κ₀ rtol = 1e-3

        @test CGL.G_real(μ, real(γ), imag(γ), κ, ξ₁, λ) ≈ [0, 0, 0, 0] atol = 1e-10
    end
end
