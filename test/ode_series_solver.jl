@testset "ode_series_solver" begin
    @testset "y'' + y = 0" begin
        # For equation y'' + y = 0
        function f(u0, t0, p; degree = 5)
            u = ArbSeries(u0[1]; degree)
            for n = 0:degree-2
                u[n+2] = -u[n] / ((n + 2) * (n + 1))
            end
            return [u]
        end

        function f(u0, p; degree = 5)
            y1 = ArbSeries(u0[1]; degree)
            y2 = ArbSeries(u0[2]; degree)

            for n = 0:degree-1
                y1[n+1] = y2[n] / (n + 1)
                y2[n+1] = -y1[n] / (n + 1)
            end

            return [y1, y2]
        end

        prob1 = ODESeriesSecondOrderProblem(
            f,
            NTuple{2,Arb}[(1, 0)],
            (Arb(0), Arb(30)),
            nothing,
        )
        prob2 = ODESeriesAutonomusProblem(f, Arb[1, 0], (Arb(0), Arb(30)), nothing)

        sol1 = ode_series_solver(prob1)
        sol2 = ode_series_solver(prob2)

        @test cos.(sol1.t) ≈ sol1[:, 1, 1] rtol = 1e-5
        @test -sin.(sol1.t) ≈ sol1[:, 1, 2] rtol = 1e-5
        # These should work once we implement a rigorous remainder term
        @test all(Arblib.overlaps.(cos.(sol1.t), sol1[:, 1, 1])) broken = true
        @test all(Arblib.overlaps.(-sin.(sol1.t), sol1[:, 1, 2])) broken = true

        @test all(Arblib.overlaps.(cos.(sol2.t), sol2[:, 1]))
        @test all(Arblib.overlaps.(-sin.(sol2.t), sol2[:, 2]))
    end
end
