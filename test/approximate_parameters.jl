@testset "approximate_parameters" begin
    paramss = [gl_params(1, 1.0, 2.3, 0.0, 0.0), gl_params(3, 1.0, 1.0, 0.0, 0.0)]
    κs = [0.49323, 0.91742]
    μs = [0.78308, 1.88590]

    ξ₁ = 30.0

    @testset "Parameters index $param_idx" for param_idx in eachindex(paramss, κs, μs)
        params, κ₀, μ₀ = paramss[param_idx], κs[param_idx], μs[param_idx]

        κ, μ = GinzburgLandauSelfSimilarSingular.approximate_parameters(κ₀, μ₀, params, ξ₁)

        # Some of the approximations are pretty bad, so we allow for a
        # large error
        @test κ ≈ κ₀ rtol = 1e-4
        @test μ ≈ μ₀ rtol = 1e-4
    end
end
