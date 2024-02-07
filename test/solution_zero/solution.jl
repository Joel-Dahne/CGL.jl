@testset "solution_zero" begin
    params = [CGL.sverak_params(Arb, 1, d) for d in [1, 3]]

    @testset "Parameters $i" for (i, (μ, γ, κ, ϵ, ξ₁, λ)) in enumerate(params)
        μ = add_error(μ, Mag(1e-8))
        κ = add_error(κ, Mag(1e-8))

        # Use a lower ξ₁. Otherwise the numerical errors are larger
        # than the enclosures.
        ξ₁ = Arb(10)

        res_capd = CGL.solution_zero_capd(μ, κ, ϵ, ξ₁, λ)
        res_float = CGL.solution_zero_float(μ, κ, ϵ, ξ₁, λ)

        res_capd_J_κ = CGL.solution_zero_jacobian_kappa_capd(μ, κ, ϵ, ξ₁, λ)
        res_float_J_κ = CGL.solution_zero_jacobian_kappa_float(μ, κ, ϵ, ξ₁, λ)

        res_capd_J_ϵ = CGL.solution_zero_jacobian_epsilon_capd(μ, κ, ϵ, ξ₁, λ)
        res_float_J_ϵ = CGL.solution_zero_jacobian_epsilon_float(μ, κ, ϵ, ξ₁, λ)

        # Check that the versions overlap
        @test all(Arblib.overlaps.(res_capd, res_float))
        @test all(Arblib.overlaps.(res_capd_J_κ, res_float_J_κ))
        @test all(Arblib.overlaps.(res_capd_J_ϵ, res_float_J_ϵ))
    end
end
