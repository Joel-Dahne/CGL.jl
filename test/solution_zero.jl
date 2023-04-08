@testset "solution_zero" begin
    paramss = [
        gl_params(1, 1.0, 2.3, 0.0, 0.0),
        gl_params(3, 1.41727, 2.3, 0.0, 0.0),
        gl_params(3, 1.41727, 2.3, 0.01, 0.02),
    ]
    κs = [0.49323, 0.45535, 0.45535]
    μs = [0.78308, 1.0, 1.0]

    ξspan = (0.0, 30.0)

    @testset "Parameters index $param_idx" for param_idx in eachindex(paramss, κs, μs)
        # Set parameters used for testing
        params, κ, μ = paramss[param_idx], κs[param_idx], μs[param_idx]
        params_arb = gl_params(Arb, params)
        u0 = SVector(μ, 0, 0, 0)

        # Numerically solve equation to have something to compare to
        prob = ODEProblem(gl_equation_real, u0, ξspan, (params, κ))
        sol = solve(prob, abstol = 1e-9, reltol = 1e-9)

        @testset "gl_taylor_expansion_real" begin
            Δξ = 0.1

            # Compute expansion at ξ0 and compare to numerical solution at
            # ξ1
            ξ0 = 0.0
            u0 = NTuple{2,Arb}[(sol(ξ0)[1], sol(ξ0)[3]), (sol(ξ0)[2], sol(ξ0)[4])]
            a_expansion, b_expansion =
                gl_taylor_expansion_real(u0, Arb(ξ0), (params_arb, Arb(κ)))

            # The tolerances are set so that it works for all choices
            # of parameters

            @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-5
            @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-4
            @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-3
            @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-3

            ξ0 = 5.0
            u0 = NTuple{2,Arb}[(sol(ξ0)[1], sol(ξ0)[3]), (sol(ξ0)[2], sol(ξ0)[4])]
            a_expansion, b_expansion =
                gl_taylor_expansion_real(u0, Arb(ξ0), (params_arb, Arb(κ)))

            @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-5
            @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-5
            @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-5
            @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-5
        end

        @testset "gl_taylor_expansion_real_autonomus vs simple" begin
            ξ0 = 0.0
            u0 = Arb[ξ0; sol(ξ0)]
            u1 = gl_taylor_expansion_real_autonomus(u0, (params_arb, Arb(κ)))
            u2 = gl_taylor_expansion_real_autonomus_simple(u0, (params_arb, Arb(κ)))
            @test all(Arblib.overlaps.(u1, u2))

            ξ0 = 5.0
            u0 = Arb[ξ0; sol(ξ0)]
            u1 = gl_taylor_expansion_real_autonomus(u0, (params_arb, Arb(κ)))
            u2 = gl_taylor_expansion_real_autonomus_simple(u0, (params_arb, Arb(κ)))

            @test all(Arblib.overlaps.(u1, u2))
        end
    end
end
