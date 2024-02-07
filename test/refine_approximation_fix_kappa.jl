@testset "refine_approximation_fix_kappa" begin
    params = [
        (0.78308, 0.49323, 0.0, CGLParams(1, 1.0, 2.3, 0.0)),
        (1.88590, 0.91742, 0.0, CGLParams(3, 1.0, 1.0, 0.0)),
    ]

    ξ₁ = 30.0

    @testset "Parameters $i" for (i, (μ₀, κ, ϵ₀, λ)) in enumerate(params)
        μ, γ, ϵ = CGL.refine_approximation_fix_kappa(μ₀, κ, ϵ₀, ξ₁, λ)

        @test μ ≈ μ₀ rtol = 1e-3
        @test ϵ ≈ ϵ₀ atol = 1e-3

        @test CGL.G(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ) ≈ [0, 0, 0, 0] atol = 1e-10
    end
end
