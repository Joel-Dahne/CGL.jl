@testset "solution_zero" begin
    paramss = [gl_params(Arb, 1, 1.0, 2.3, 0.0, 0.0), gl_params(Arb, 3, 1.0, 1.0, 0.0, 0.0)]
    κs = add_error.(Arb[0.49323, 0.91742], Mag(1e-9))
    μs = add_error.(Arb[0.78308, 1.88590], Mag(1e-9))

    ξ₁ = Arb(10.0)

    @testset "Parameters index $param_idx" for param_idx in eachindex(paramss, κs, μs)
        params, κ, μ = paramss[param_idx], κs[param_idx], μs[param_idx]

        res_capd = GinzburgLandauSelfSimilarSingular.solution_zero_capd(μ, κ, ξ₁, params)
        res_float = GinzburgLandauSelfSimilarSingular.solution_zero_float(μ, κ, ξ₁, params)

        @test all(Arblib.overlaps.(res_capd, res_float))
    end
end
