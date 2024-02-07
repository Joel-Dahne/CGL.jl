@testset "G" begin
    G = CGL.G
    G_jacobian_kappa = CGL.G_jacobian_kappa
    G_jacobian_epsilon = CGL.G_jacobian_epsilon

    params = [CGL.sverak_params(Float64, 1, d) for d in [1, 3]]

    @testset "Parameters $i" for (i, (μ, γ, κ, ϵ, ξ₁, λ)) in enumerate(params)
        J_kappa = G_jacobian_kappa(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)
        J_kappa_Arb = G_jacobian_kappa(
            Arb(μ),
            Arb(real(γ)),
            Arb(imag(γ)),
            Arb(κ),
            Arb(ϵ),
            Arb(ξ₁),
            CGLParams{Arb}(λ),
        )

        J_epsilon = G_jacobian_epsilon(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ)
        J_epsilon_Arb = G_jacobian_epsilon(
            Arb(μ),
            Arb(real(γ)),
            Arb(imag(γ)),
            Arb(κ),
            Arb(ϵ),
            Arb(ξ₁),
            CGLParams{Arb}(λ),
        )

        @test J_kappa ≈ Float64.(J_kappa_Arb) rtol = 1e-4
        @test J_epsilon ≈ Float64.(J_epsilon_Arb) rtol = 1e-4

        # Derivative w.r.t. to μ should overlap
        @test all(Arblib.overlaps.(J_kappa_Arb[:, 1], J_epsilon_Arb[:, 1]))
        # Derivative w.r.t. γ should be identical
        @test isequal(J_kappa_Arb[:, 2:3], J_epsilon_Arb[:, 2:3])

        # Function for computing derivative using finite differences.
        fdm = central_fdm(5, 1; factor = 1e5)

        @test fdm(μ -> G(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ), μ) ≈ J_kappa[:, 1] rtol = 1e-4
        @test fdm(rγ -> G(μ, rγ, imag(γ), κ, ϵ, ξ₁, λ), real(γ)) ≈ J_kappa[:, 2] rtol = 1e-4
        @test fdm(iγ -> G(μ, real(γ), iγ, κ, ϵ, ξ₁, λ), imag(γ)) ≈ J_kappa[:, 3] rtol = 1e-4

        @test fdm(κ -> G(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ), κ) ≈ J_kappa[:, 4] rtol = 1e-2
        @test fdm(ϵ -> G(μ, real(γ), imag(γ), κ, ϵ, ξ₁, λ), ϵ) ≈ J_epsilon[:, 4] rtol = 1e-2
    end

    @testset "Isolate zero" begin
        @testset "Parameters $i" for (i, (μ, γ, κ, ϵ, ξ₁, λ)) in enumerate(params)
            # Ball in which we want to verify root
            x₀_F64 = SVector(μ, real(γ), imag(γ), κ)
            x₀ = Arb.(x₀_F64)
            x = add_error.(x₀, Mag(1e-6))

            # Verify root
            xx = let ϵ = Arb(ϵ), ξ₁ = Arb(ξ₁), λ = CGLParams{Arb}(λ)
                G_x = x -> G(x..., ϵ, ξ₁, λ)
                dG_x = x -> G_jacobian_kappa(x..., ϵ, ξ₁, λ)

                CGL.verify_and_refine_root(G_x, dG_x, x, max_iterations = 5)
            end

            @test all(isfinite, xx)
        end
    end
end
