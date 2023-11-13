@testset "G_solve" begin
    @testset "Parameters (i, d)" for (i, d) in [(1, 1), (2, 1)]
        param1 = CGL._params(Arb, i, d, ξ₁ = Arb(30))
        param2 = CGL._params(Arb, i, d, ξ₁ = Arb(35))

        res1 = CGL.G_solve(param1...)
        res2 = CGL.G_solve(param2...)

        @test all(isfinite, res1)
        @test all(isfinite, res2)
        @test Arblib.overlaps(res1[1], res2[1])
        @test Arblib.overlaps(res1[4], res2[4])
    end
end
