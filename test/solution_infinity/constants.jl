@testset "solution_infinity_constants" begin
    params = [
        (Arb(0.493223), gl_params(Arb, 1, 1.0, 2.3, 0.0, 0.0)),
        (Arb(0.917383), gl_params(Arb, 3, 1.0, 1.0, 0.0, 0.0)),
        (Arb(0.917383), gl_params(Arb, 3, 1.0, 1.0, 0.01, 0.02)),
    ]

    ξ₁ = Arb(10)

    @testset "C_T1" begin
        # We don't have an obvious way to check this inequality. For
        # now we simply check that C₁ is finite and C₂ it is
        # increasing in v.
        for (κ, λ) in params
            C1₁, C1₂ = GinzburgLandauSelfSimilarSingular.C_T1(Arb(0.1), κ, λ, ξ₁)
            C2₁, C2₂ = GinzburgLandauSelfSimilarSingular.C_T1(Arb(0.2), κ, λ, ξ₁)

            @test isfinite(C1₁)
            @test isfinite(C2₁)
            @test C1₂ < C2₂
        end
    end
end
