@testset "Q_infinity" begin
    params = [CGL.sverak_params.(Arb, 1, d) for d in [1, 3]]

    @testset "Q_infinity" begin
        @testset "Parameters $i" for (i, (μ, γ, κ, ϵ, ξ₁, λ)) in enumerate(params)
            res_F64 = CGL.Q_infinity(
                Complex{Float64}(γ),
                Float64(κ),
                Float64(ϵ),
                Float64(ξ₁),
                CGLParams{Float64}(λ),
            )
            res_Arb = CGL.Q_infinity(γ, κ, ϵ, ξ₁, λ)

            @test res_F64 ≈ ComplexF64.(res_Arb) rtol = 1e-5
        end
    end

    @testset "Q_infinity_jacobian_kappa" begin
        @testset "Parameters $i" for (i, (μ, γ, κ, ϵ, ξ₁, λ)) in enumerate(params)
            res_J_F64 = CGL.Q_infinity_jacobian_kappa(
                Complex{Float64}(γ),
                Float64(κ),
                Float64(ϵ),
                Float64(ξ₁),
                CGLParams{Float64}(λ),
            )

            res_J_Arb = CGL.Q_infinity_jacobian_kappa(γ, κ, ϵ, ξ₁, λ)

            # For the first column we get good enclosures
            @test res_J_F64[:, 1] ≈ ComplexF64.(res_J_Arb[:, 1]) rtol = 1e-3

            # For the second column we get very bad enclosures and
            # hence not very good agreement. We only check that the
            # F64 approximation is inside the compute enclosure.
            @test all(Arblib.contains.(res_J_Arb[:, 2], Acb.(res_J_F64[:, 2])))
        end
    end

    @testset "Q_infinity_jacobian_epsilon" begin
        @testset "Parameters $i" for (i, (μ, γ, κ, ϵ, ξ₁, λ)) in enumerate(params)
            res_J_F64 = CGL.Q_infinity_jacobian_epsilon(
                Complex{Float64}(γ),
                Float64(κ),
                Float64(ϵ),
                Float64(ξ₁),
                CGLParams{Float64}(λ),
            )

            res_J_Arb = CGL.Q_infinity_jacobian_epsilon(γ, κ, ϵ, ξ₁, λ)

            # For the first column we get good enclosures
            @test res_J_F64[:, 1] ≈ ComplexF64.(res_J_Arb[:, 1]) rtol = 1e-3

            # For the second column we get very bad enclosures and
            # hence not very good agreement. We only check that the
            # F64 approximation is inside the compute enclosure.
            @test all(Arblib.contains.(res_J_Arb[:, 2], Acb.(res_J_F64[:, 2])))
        end
    end
end
