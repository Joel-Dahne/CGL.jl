@testset "G_solve" begin
    @testset "Parameters (i, d)" for (i, d) in [(1, 1), (2, 1)]
        μ1, γ1, κ1, ϵ1, ξ₁1, λ1 = CGL.sverak_params(Arb, i, d, ξ₁ = Arb(30))
        μ2, γ2, κ2, ϵ2, ξ₁2, λ2 = CGL.sverak_params(Arb, i, d, ξ₁ = Arb(35))

        res1 = CGL.G_solve_fix_epsilon_alt(μ1, real(γ1), imag(γ1), κ1, ϵ1, ξ₁1, λ1)
        res2 = CGL.G_solve_fix_epsilon_alt(μ2, real(γ2), imag(γ2), κ2, ϵ2, ξ₁2, λ2)

        @test all(isfinite, res1)
        @test all(isfinite, res2)
        @test Arblib.overlaps(res1[1], res2[1])
        @test Arblib.overlaps(res1[4], res2[4])
    end
end
