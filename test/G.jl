@testset "G" begin
    G = GinzburgLandauSelfSimilarSingular.G_real
    G_jacobian = GinzburgLandauSelfSimilarSingular.G_jacobian_real

    params = [
        (0.783077, 1.79202, -1.10431, 0.493223, gl_params(1, 1.0, 2.3, 0.0, 0.0)),
        (1.88576, 1.71299, -1.49238, 0.917383, gl_params(3, 1.0, 1.0, 0.0, 0.0)),
    ]

    ξ₁ = 30.0

    @testset "Parameters $i" for (i, (μ, γ_real, γ_imag, κ, λ)) in enumerate(params)
        res = G(μ, γ_real, γ_imag, κ, ξ₁, λ)
        res_jacobian = G_jacobian(μ, γ_real, γ_imag, κ, ξ₁, λ)

        # Compare central difference with jacobian. It is highly
        # oscillating in κ so we have a high tolerance for that
        # derivative.
        Δ = 0.0001
        @test (G(μ + Δ, γ_real, γ_imag, κ, ξ₁, λ) - G(μ - Δ, γ_real, γ_imag, κ, ξ₁, λ)) /
              2Δ ≈ res_jacobian[:, 1] atol = Δ
        @test (G(μ, γ_real + Δ, γ_imag, κ, ξ₁, λ) - G(μ, γ_real - Δ, γ_imag, κ, ξ₁, λ)) /
              2Δ ≈ res_jacobian[:, 2] atol = Δ
        @test (G(μ, γ_real, γ_imag + Δ, κ, ξ₁, λ) - G(μ, γ_real, γ_imag - Δ, κ, ξ₁, λ)) /
              2Δ ≈ res_jacobian[:, 3] atol = Δ
        @test (G(μ, γ_real, γ_imag, κ + Δ, ξ₁, λ) - G(μ, γ_real, γ_imag, κ - Δ, ξ₁, λ)) /
              2Δ ≈ res_jacobian[:, 4] atol = 0.05
    end

    @testset "Isolate zero" begin
        @testset "Parameters $i" for (i, (μ, γ_real, γ_imag, κ, λ)) in enumerate(params)
            # Find approximate zero
            μ₀, γ₀, κ₀ =
                GinzburgLandauSelfSimilarSingular.approximate_parameters(μ, κ, ξ₁, λ)

            @test μ ≈ μ₀ rtol = 1e-3
            @test γ_real ≈ real(γ₀) rtol = 1e-3
            @test γ_imag ≈ imag(γ₀) rtol = 1e-3
            @test κ ≈ κ₀ rtol = 1e-3

            # Ball in which we want to verify root
            x₀_F64 = SVector(μ₀, real(γ₀), imag(γ₀), κ₀)
            x₀ = Arb.(x₀_F64)
            x = add_error.(x₀, Mag(1e-6))

            # Verify root
            xx = let ξ₁ = Arb(ξ₁), λ = gl_params(Arb, λ)
                G_x = x -> G(x[1], x[2], x[3], x[4], ξ₁, λ)
                dG_x = x -> G_jacobian(x[1], x[2], x[3], x[4], ξ₁, λ)

                GinzburgLandauSelfSimilarSingular.verify_and_refine_root(
                    G_x,
                    dG_x,
                    x,
                    max_iterations = 5,
                )
            end

            @test all(isfinite, xx) broken = λ.d == 3
        end
    end
end
