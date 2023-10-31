@testset "G" begin
    G = GinzburgLandauSelfSimilarSingular.G_real
    G_jacobian = GinzburgLandauSelfSimilarSingular.G_jacobian_real

    params = [GinzburgLandauSelfSimilarSingular._params.(Float64, 1, d) for d in [1, 3]]

    @testset "Parameters $i" for (i, (μ, γ, κ, ξ₁, λ)) in enumerate(params)
        res = G_jacobian(μ, real(γ), imag(γ), κ, ξ₁, λ)

        # Function for computing derivative using finite differences.
        fdm = central_fdm(5, 1; factor = 1e5)

        @test fdm(μ -> G(μ, real(γ), imag(γ), κ, ξ₁, λ), μ) ≈ res[:, 1] atol = 1e-5
        @test fdm(rγ -> G(μ, rγ, imag(γ), κ, ξ₁, λ), real(γ)) ≈ res[:, 2] atol = 1e-5
        @test fdm(iγ -> G(μ, real(γ), iγ, κ, ξ₁, λ), imag(γ)) ≈ res[:, 3] atol = 1e-5
        @test fdm(κ -> G(μ, real(γ), imag(γ), κ, ξ₁, λ), κ) ≈ res[:, 4] atol = 1e-4
    end

    @testset "Isolate zero" begin
        @testset "Parameters $i" for (i, (μ, γ, κ, ξ₁, λ)) in enumerate(params)
            # Ball in which we want to verify root
            x₀_F64 = SVector(μ, real(γ), imag(γ), κ)
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
