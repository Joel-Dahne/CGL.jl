@testset "check_existence" begin
    paramss = [
        gl_params(1, Arb(1.0), 2.3, 0.0, 0.0),
        gl_params(3, Arb(1.41727), 1.0, 0.0, 0.0),
        gl_params(3, Arb(1.41727), 1.0, 0.01, 0.02),
    ]

    κs = add_error.(Arb[0.49323, 0.45535, 0.45535], Mag(1e-10))
    μs = add_error.(Arb[0.78308, 1.0, 1.0], Mag(1e-10))

    ξ₁ = Arb(10)
    v = Arb(0.1)

    @testset "Parameters index $param_idx" for param_idx in eachindex(paramss, κs, μs)
        # Set parameters used for testing
        params, κ, μ = paramss[param_idx], κs[param_idx], μs[param_idx]

        if params.d == 1
            sucess1, γ1 = GinzburgLandauSelfSimilarSingular.check_existence_fixed_point(
                κ,
                μ,
                params,
                ξ₁,
                v,
                non_rigorous = false,
            )

            sucess2, γ2 = GinzburgLandauSelfSimilarSingular.check_existence_fixed_point(
                κ,
                μ,
                params,
                ξ₁,
                v,
                non_rigorous = true,
            )

            @test sucess1
            @test sucess2
            @test Arblib.overlaps(γ1, γ2)
        else
            @test_throws ArgumentError GinzburgLandauSelfSimilarSingular.check_existence_fixed_point(
                κ,
                μ,
                params,
                ξ₁,
                v,
            )

            sucess2, γ2 = GinzburgLandauSelfSimilarSingular.check_existence_fixed_point(
                κ,
                μ,
                params,
                ξ₁,
                v,
                non_rigorous = true,
            )

            @test sucess2
            @test isfinite(γ2)
        end
    end
end
