@testset "G" begin
    params = [
        (0.783077, 1.79202, -1.10431, 0.493223, gl_params(1, 1.0, 2.3, 0.0, 0.0)),
        (1.88576, 1.71299, -1.49238, 0.917383, gl_params(3, 1.0, 1.0, 0.0, 0.0)),
    ]

    ξ₁ = 40.0

    @testset "Parameters $i" for (i, (μ, γ_real, γ_imag, κ, λ)) in enumerate(params)
        G = GinzburgLandauSelfSimilarSingular.G_real

        res = G(μ, γ_real, γ_imag, κ, ξ₁, λ)
        res_jacobian =
            GinzburgLandauSelfSimilarSingular.G_jacobian_real(μ, γ_real, γ_imag, κ, ξ₁, λ)

        # Check that we get a decent first order approximation
        Δ = 0.01
        @test_broken G(μ + Δ, γ_real, γ_imag, κ, ξ₁, λ) ≈ res + Δ * res_jacobian[:, 1]
        @test G(μ, γ_real + Δ, γ_imag, κ, ξ₁, λ) ≈ res + Δ * res_jacobian[:, 2]
        @test G(μ, γ_real, γ_imag + Δ, κ, ξ₁, λ) ≈ res + Δ * res_jacobian[:, 3]
        @test_broken G(μ, γ_real, γ_imag, κ + Δ, ξ₁, λ) ≈ res + Δ * res_jacobian[:, 4]
    end
end
