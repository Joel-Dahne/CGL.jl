@testset "solution_zero_equation" begin
    params = [
        (0.783077, 0.493223, gl_params(1, 1.0, 2.3, 0.0, 0.0)),
        (1.88576, 0.917383, gl_params(3, 1.0, 1.0, 0.0, 0.0)),
        (1.88576, 0.917383, gl_params(3, 1.0, 1.0, 0.01, 0.02)),
    ]

    ξspan = (0.0, 30.0)

    @testset "Parameters $i" for (i, (μ, κ, λ)) in enumerate(params)
        λ_Arb = gl_params(Arb, λ)
        u0 = SVector(μ, 0, 0, 0)

        # Numerically solve equation to have something to compare to
        prob = ODEProblem(gl_equation_real_system_ode, u0, ξspan, (κ, λ))
        sol = solve(prob, abstol = 1e-9, reltol = 1e-9)

        @testset "gl_equation_real_taylor_expansion" begin
            Δξ = 0.1

            # Compute expansion at ξ0 and compare to numerical solution at
            # ξ1
            ξ0 = 0.0
            u0 = NTuple{2,Arb}[(sol(ξ0)[1], sol(ξ0)[3]), (sol(ξ0)[2], sol(ξ0)[4])]
            a_expansion, b_expansion =
                gl_equation_real_taylor_expansion(u0, Arb(κ), Arb(ξ0), λ_Arb)

            # The tolerances are set so that it works for all choices
            # of parameters

            @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-5
            @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-4
            @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-3
            @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-3

            ξ0 = 5.0
            u0 = NTuple{2,Arb}[(sol(ξ0)[1], sol(ξ0)[3]), (sol(ξ0)[2], sol(ξ0)[4])]
            a_expansion, b_expansion =
                gl_equation_real_taylor_expansion(u0, Arb(κ), Arb(ξ0), λ_Arb)

            @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-5
            @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-5
            @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-4
            @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-5
        end

        @testset "gl_equation_real_system_autonomus_taylor_expansion vs simple" begin
            ξ0 = 0.0
            u0 = Arb[sol(ξ0); ξ0]
            u1 = gl_equation_real_system_autonomus_taylor_expansion(u0, Arb(κ), λ_Arb)
            u2 =
                gl_equation_real_system_autonomus_taylor_expansion_simple(u0, Arb(κ), λ_Arb)
            @test isfinite.(u1) == isfinite.(u2)
            @test all(Arblib.overlaps.(u1, u2))

            # Check that the radius is not too different for them
            r1 = [radius(u1[i][j]) for i in eachindex(u1), j = 0:Arblib.degree(u1[1])]
            r2 = [radius(u2[i][j]) for i in eachindex(u2), j = 0:Arblib.degree(u2[1])]
            @test all(r -> 1 / 8 < Float64(r) < 8, filter(isfinite, r1 ./ r2))

            ξ0 = 5.0
            u0 = Arb[sol(ξ0); ξ0]
            u1 = gl_equation_real_system_autonomus_taylor_expansion(u0, Arb(κ), λ_Arb)
            u2 =
                gl_equation_real_system_autonomus_taylor_expansion_simple(u0, Arb(κ), λ_Arb)

            @test isfinite.(u1) == isfinite.(u2)
            @test all(Arblib.overlaps.(u1, u2))

            # Check that the radius is not too different for them
            r1 = [radius(u1[i][j]) for i in eachindex(u1), j = 0:Arblib.degree(u1[1])]
            r2 = [radius(u2[i][j]) for i in eachindex(u2), j = 0:Arblib.degree(u2[1])]

            @test all(r -> 1 / 16 < Float64(r) < 16, filter(isfinite, r1 ./ r2))
        end
    end
end
