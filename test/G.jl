@testset "G" begin
    G = CGL.G_real
    G_jacobian = CGL.G_jacobian_real
    G_jacobian_epsilon = CGL.G_jacobian_epsilon_real

    params = [CGL.sverak_params(Float64, 1, d) for d in [1, 3]]

    @testset "Parameters $i" for (i, (μ, γ, κ, ξ₁, λ)) in enumerate(params)
        res = G_jacobian(μ, real(γ), imag(γ), κ, ξ₁, λ)
        res_Arb = G_jacobian(
            Arb(μ),
            Arb(real(γ)),
            Arb(imag(γ)),
            Arb(κ),
            Arb(ξ₁),
            CGLParams{Arb}(λ),
        )

        @test res ≈ Float64.(res_Arb) rtol = 1e-4

        # Function for computing derivative using finite differences.
        fdm = central_fdm(5, 1; factor = 1e5)

        @test fdm(μ -> G(μ, real(γ), imag(γ), κ, ξ₁, λ), μ) ≈ res[:, 1] rtol = 1e-4
        @test fdm(rγ -> G(μ, rγ, imag(γ), κ, ξ₁, λ), real(γ)) ≈ res[:, 2] rtol = 1e-4
        @test fdm(iγ -> G(μ, real(γ), iγ, κ, ξ₁, λ), imag(γ)) ≈ res[:, 3] rtol = 1e-4
        @test fdm(κ -> G(μ, real(γ), imag(γ), κ, ξ₁, λ), κ) ≈ res[:, 4] rtol = 1e-2

        if λ.d == 1
            res_epsilon_Arb = G_jacobian_epsilon(
                Arb(μ),
                Arb(real(γ)),
                Arb(imag(γ)),
                Arb(κ),
                Arb(ξ₁),
                CGLParams{Arb}(λ),
            )
            res_epsilon = Float64.(res_epsilon_Arb) # No Float64 version

            # Derivative w.r.t. to μ should overlap
            @test all(Arblib.overlaps.(res_Arb[:, 1], res_epsilon_Arb[:, 1]))
            # Derivative w.r.t. γ should be identical
            @test isequal(res_Arb[:, 2:3], res_epsilon_Arb[:, 2:3])

            @test fdm(ϵ -> G(μ, real(γ), imag(γ), κ, ξ₁, CGLParams(λ; ϵ)), λ.ϵ) ≈
                  res_epsilon[:, 4] rtol = 1e-4
        end
    end

    @testset "Isolate zero" begin
        @testset "Parameters $i" for (i, (μ, γ, κ, ξ₁, λ)) in enumerate(params)
            # Ball in which we want to verify root
            x₀_F64 = SVector(μ, real(γ), imag(γ), κ)
            x₀ = Arb.(x₀_F64)
            x = add_error.(x₀, Mag(1e-6))

            # Verify root
            xx = let ξ₁ = Arb(ξ₁), λ = CGLParams{Arb}(λ)
                G_x = x -> G(x..., ξ₁, λ)
                dG_x = x -> G_jacobian(x..., ξ₁, λ)

                CGL.verify_and_refine_root(G_x, dG_x, x, max_iterations = 5)
            end

            @test all(isfinite, xx)
        end
    end
end
