@testset "solution_zero" begin
    params = [GinzburgLandauSelfSimilarSingular._params.(Arb, 1, d) for d in [1, 3]]

    @testset "Parameters $i" for (i, (μ, γ, κ, ξ₁, λ)) in enumerate(params)
        μ = add_error(μ, Mag(1e-8))
        κ = add_error(κ, Mag(1e-8))

        # Use a lower ξ₁. Otherwise the numerical errors are larger
        # than the enclosures.
        ξ₁ = Arb(10)

        res_capd = GinzburgLandauSelfSimilarSingular.solution_zero_capd(μ, κ, ξ₁, λ)
        res_float = GinzburgLandauSelfSimilarSingular.solution_zero_float(μ, κ, ξ₁, λ)

        res_capd_2, res_capd_jac =
            GinzburgLandauSelfSimilarSingular.solution_zero_jacobian_capd(μ, κ, ξ₁, λ)
        res_float_2, res_float_jac =
            GinzburgLandauSelfSimilarSingular.solution_zero_jacobian_float(μ, κ, ξ₁, λ)

        # Check that all versions overlap
        @test all(Arblib.overlaps.(res_capd, res_float))
        @test all(Arblib.overlaps.(res_capd, res_capd_2))
        @test all(Arblib.overlaps.(res_capd, res_float_2))
        @test all(Arblib.overlaps.(res_capd_2, res_float))
        @test all(Arblib.overlaps.(res_capd_2, res_float_2))

        @test all(Arblib.overlaps.(res_capd_jac, res_float_jac)) broken = (λ.d != 1)
        if λ.d != 1
            @test res_float_jac \ res_capd_jac ≈ LinearAlgebra.I atol = 1e-2
        end
    end
end
