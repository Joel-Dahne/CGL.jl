@testset "Q_zero_equation" begin
    params = [
        (0.783077, 0.493223, 0.0, CGLParams(1, 1.0, 2.3, 0.0)),
        (0.783077, 0.493223, 0.1, CGLParams(1, 1.0, 2.3, 0.01)),
        (1.88576, 0.917383, 0.0, CGLParams(3, 1.0, 1.0, 0.0)),
        (1.88576, 0.917383, 0.1, CGLParams(3, 1.0, 1.0, 0.01)),
    ]

    @testset "Parameters $i" for (i, (μ, κ, ϵ, λ)) in enumerate(params)
        λ_Arb = CGLParams{Arb}(λ)
        u0 = SVector(μ, 0, 0, 0)

        # Numerically solve equation to have something to compare to
        sol = CGL.Q_zero_float_curve(μ, κ, ϵ, 10.0, λ)

        @testset "cgl_equation_real_taylor" begin
            Δξ = 0.1

            # Compute expansion at ξ0 and compare to numerical solution at
            # ξ1
            ξ0 = 0.0
            u0 = NTuple{2,Arb}[(sol(ξ0)[1], sol(ξ0)[3]), (sol(ξ0)[2], sol(ξ0)[4])]
            a_expansion, b_expansion = CGL.cgl_equation_real_taylor(
                u0,
                Arb(κ),
                Arb(ϵ),
                Arb(ξ0),
                λ_Arb,
                degree = 10,
            )

            @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-8
            @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-8
            @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-8
            @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-8

            ξ0 = 5.0
            u0 = NTuple{2,Arb}[(sol(ξ0)[1], sol(ξ0)[3]), (sol(ξ0)[2], sol(ξ0)[4])]
            a_expansion, b_expansion = CGL.cgl_equation_real_taylor(
                u0,
                Arb(κ),
                Arb(ϵ),
                Arb(ξ0),
                λ_Arb,
                degree = 10,
            )

            @test a_expansion(Δξ) ≈ sol(ξ0 + Δξ)[1] rtol = 1e-8
            @test b_expansion(Δξ) ≈ sol(ξ0 + Δξ)[2] rtol = 1e-8
            @test Arblib.derivative(a_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[3] rtol = 1e-8
            @test Arblib.derivative(b_expansion)(Δξ) ≈ sol(ξ0 + Δξ)[4] rtol = 1e-8
        end
    end
end
