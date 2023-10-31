@testset "functions" begin
    ξ = Arb(30)
    κ = Arb(0.493223)
    λ = gl_params(Arb, 1, 1.0, 2.3, 0.0, 0.0)

    @testset "P" begin
        @test Arblib.overlaps(P(ArbSeries((ξ, 1)), (λ, κ))[1], P_dξ(ξ, (λ, κ)))

        @test Arblib.overlaps(P(ξ, (λ, ArbSeries((κ, 1))))[1], P_dκ(ξ, (λ, κ)))

        @test CGL.P_asym_approx(ξ, (λ, κ)) ≈ P(ξ, (λ, κ)) rtol = 1e-20

        @test CGL.P_dξ_asym_approx(ξ, (λ, κ)) ≈ P_dξ(ξ, (λ, κ)) rtol = 1e-20
    end

    @testset "E" begin
        @test Arblib.overlaps(E(ArbSeries((ξ, 1)), (λ, κ))[1], E_dξ(ξ, (λ, κ)))

        @test Arblib.overlaps(E(ξ, (λ, ArbSeries((κ, 1))))[1], E_dκ(ξ, (λ, κ)))
    end
end
