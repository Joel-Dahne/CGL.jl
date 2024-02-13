@testset "special-functions" begin
    @testset "rising" begin
        @test rising(Arb(5), Arb(3)) == 210
        @test rising(Arb(5), 3) == 210
        @test rising(ArbSeries((5, 1)), 3) == ArbSeries((210, 107))
        @test rising(5.0, 3.0) == 210
        @test rising(5.0, 3) == 210
        @test rising(5, 3) == 210
    end
end
