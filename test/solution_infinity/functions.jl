@testset "functions" begin
    ξ = Arb(30)
    κ = Arb(0.493223)
    λ = gl_params(Arb, 1, 1.0, 2.3, 0.1, 0.2)
    _, _, c = CGL._abc(κ, λ)

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

    @testset "J_P" begin
        @test Arblib.overlaps(
            J_P(ξ, (λ, κ)),
            CGL.B_W(κ, λ) * P(ξ, (λ, κ)) * exp(-c * ξ^2) * ξ^(λ.d - 1),
        )

        η = Arb(40)
        @test Arblib.overlaps(
            E(ξ, (λ, κ)) * J_P(η, (λ, κ)),
            -(1 + im * λ.δ) * K(ξ, η, (λ, κ)),
        )
    end

    @testset "J_E" begin
        @test Arblib.overlaps(
            J_E(ξ, (λ, κ)),
            CGL.B_W(κ, λ) * E(ξ, (λ, κ)) * exp(-c * ξ^2) * ξ^(λ.d - 1),
        )

        η = Arb(20)
        @test Arblib.overlaps(
            P(ξ, (λ, κ)) * J_E(η, (λ, κ)),
            -(1 + im * λ.δ) * K(ξ, η, (λ, κ)),
        )
    end
end
