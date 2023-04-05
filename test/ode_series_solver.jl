@testset "ode_series_solver" begin
    # For equation y'' + y = 0
    function f(u0, t0, p; degree = 5)
        u = ArbSeries(u0[1]; degree)
        for n = 0:degree-2
            u[n+2] = -u[n] / ((n + 2) * (n + 1))
        end
        return [u]
    end

    prob = ODESeriesSecondOrderProblem(f, NTuple{2,Arb}[(1, 0)], (Arb(0), Arb(30)), nothing)

    sol = ode_series_solver(prob)

    u = getindex.(sol.u, 1)

    @test cos.(sol.t) ≈ getindex.(u, 1) rtol = 1e-5
    @test -sin.(sol.t) ≈ getindex.(u, 2) rtol = 1e-5
    # These should work once we implement a rigorous remainder term
    @test all(Arblib.overlaps.(cos.(sol.t), getindex.(u, 1))) broken = true
    @test all(Arblib.overlaps.(-sin.(sol.t), getindex.(u, 2))) broken = true
end
