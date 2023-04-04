@testset "solution_zero" begin
    # Set parameters used for testing
    params = gl_params(1, 1.0, 2.3, 0.0, 0.0)
    κ, μ = 0.49323, 0.78308
    ξspan = (0.0, 30.0)
    u0 = SVector(μ, 0, 0, 0)

    # Numerically solve equation to have something to compare to
    prob = ODEProblem(gl_equation_real, u0, ξspan, (params, κ))
    sol = solve(prob, abstol = 1e-9, reltol = 1e-9)

    @testset "gl_taylor_expansion_real_1" begin
        params_arb = gl_params(params.d, Arb(params.ω), params.σ, params.ϵ, params.δ)
        Δξ = 0.1

        # Compute expansion at ξ0 and compare to numerical solution at
        # ξ1
        ξ0 = 0.0
        a_expansion, b_expansion =
            gl_taylor_expansion_real(Arb.(sol(ξ0)), Arb(ξ0), (params_arb, Arb(κ)))


        @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-6
        @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-5
        @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-3
        @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-3

        ξ0 = 5.0
        a_expansion, b_expansion =
            gl_taylor_expansion_real(Arb.(sol(ξ0)), Arb(ξ0), (params_arb, Arb(κ)))

        @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-8
        @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-8
        @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-6
        @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-6
    end

    # Set parameters used for testing
    params = gl_params(3, 1.41727, 2.3, 0.0, 0.0)
    κ, μ = 0.45535, 1.0
    ξspan = (0.0, 30.0)
    u0 = SVector(μ, 0, 0, 0)

    # Numerically solve equation to have something to compare to
    prob = ODEProblem(gl_equation_real, u0, ξspan, (params, κ))
    sol = solve(prob, abstol = 1e-9, reltol = 1e-9)

    @testset "gl_taylor_expansion_real_2" begin
        params_arb = gl_params(params.d, Arb(params.ω), params.σ, params.ϵ, params.δ)
        Δξ = 0.1

        # Compute expansion at ξ0 and compare to numerical solution at
        # ξ0 + Δξ
        ξ0 = 0.0
        a_expansion, b_expansion =
            gl_taylor_expansion_real(Arb.(sol(ξ0)), Arb(ξ0), (params_arb, Arb(κ)))

        @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-6
        @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-5
        @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-3
        @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-3

        ξ0 = 5.0
        a_expansion, b_expansion =
            gl_taylor_expansion_real(Arb.(sol(ξ0)), Arb(ξ0), (params_arb, Arb(κ)))

        @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-8
        @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-8
        @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-5
        @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-5
    end

    # Set parameters used for testing
    params = gl_params(3, 1.41727, 2.3, 0.01, 0.02)
    κ, μ = 0.45535, 1.0
    ξspan = (0.0, 30.0)
    u0 = SVector(μ, 0, 0, 0)

    # Numerically solve equation to have something to compare to
    prob = ODEProblem(gl_equation_real, u0, ξspan, (params, κ))
    sol = solve(prob, abstol = 1e-9, reltol = 1e-9)

    @testset "gl_taylor_expansion_real_3" begin
        params_arb = gl_params(params.d, Arb(params.ω), params.σ, params.ϵ, params.δ)
        Δξ = 0.1

        # Compute expansion at ξ0 and compare to numerical solution at
        # ξ1
        ξ0 = 0.0
        a_expansion, b_expansion =
            gl_taylor_expansion_real(Arb.(sol(ξ0)), Arb(ξ0), (params_arb, Arb(κ)))


        @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-6
        @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-4
        @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-4
        @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-4

        ξ0 = 5.0
        a_expansion, b_expansion =
            gl_taylor_expansion_real(Arb.(sol(ξ0)), Arb(ξ0), (params_arb, Arb(κ)))

        @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-8
        @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-8
        @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-5
        @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-5
    end
end
