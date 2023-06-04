@testset "solution_zero" begin
    paramss = [
        gl_params(1, 1.0, 2.3, 0.0, 0.0),
        gl_params(3, 1.41727, 2.3, 0.0, 0.0),
        gl_params(3, 1.41727, 2.3, 0.01, 0.02),
    ]
    κs = [0.49323, 0.45535, 0.45535]
    μs = [0.78308, 1.0, 1.0]

    ξ₁ = 30.0

    @testset "Parameters index $param_idx" for param_idx in eachindex(paramss, κs, μs)
        # Set parameters used for testing
        params, κ, μ = paramss[param_idx], κs[param_idx], μs[param_idx]

        if params.d == 1
            res_capd =
                GinzburgLandauSelfSimilarSingular.solution_zero_capd(κ, μ, params, ξ₁)
            res_float =
                GinzburgLandauSelfSimilarSingular.solution_zero_float(κ, μ, params, ξ₁)

            @test res_capd ≈ res_float
        else
            # The CAPD program throws an exception in this case
            @test_throws ArgumentError GinzburgLandauSelfSimilarSingular.solution_zero_capd(
                κ,
                μ,
                params,
                ξ₁,
            )
        end
    end
end
