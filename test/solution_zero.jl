@testset "solution_zero" begin
    paramss = [
        gl_params(1, Arb(1.0), 2.3, 0.0, 0.0),
        gl_params(3, Arb(1.41727), 2.3, 0.0, 0.0),
        gl_params(3, Arb(1.41727), 2.3, 0.01, 0.02),
    ]
    κs = add_error.(Arb[0.49323, 0.45535, 0.45535], Mag(1e-10))
    μs = add_error.(Arb[0.78308, 1.0, 1.0], Mag(1e-10))

    ξ₁ = Arb(10.0)

    @testset "Parameters index $param_idx" for param_idx in eachindex(paramss, κs, μs)
        params, κ, μ = paramss[param_idx], κs[param_idx], μs[param_idx]

        if params.d == 1
            res_capd =
                GinzburgLandauSelfSimilarSingular.solution_zero_capd(κ, μ, params, ξ₁)
            res_float =
                GinzburgLandauSelfSimilarSingular.solution_zero_float(κ, μ, params, ξ₁)

            @test all(Arblib.overlaps.(res_capd, res_float))
        else
            # The CAPD program throws an exception in this case
            @test_throws ArgumentError GinzburgLandauSelfSimilarSingular.solution_zero_capd(
                κ,
                μ,
                params,
                ξ₁,
            )

            res_float =
                GinzburgLandauSelfSimilarSingular.solution_zero_float(κ, μ, params, ξ₁)

            @test all(isfinite, res_float)
        end
    end
end
