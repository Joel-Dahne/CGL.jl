@testset "solution_infinity_constants" begin
    params = [
        (Arb(0.493223), gl_params(Arb, 1, 1.0, 2.3, 0.0, 0.0)),
        (Arb(0.917383), gl_params(Arb, 3, 1.0, 1.0, 0.0, 0.0)),
        (Arb(0.917383), gl_params(Arb, 3, 1.0, 1.0, 0.01, 0.02)),
    ]

    ξ₁ = Arb(10)

    @testset "C_P" begin
        for (κ, λ) in params
            C = GinzburgLandauSelfSimilarSingular.C_P(κ, λ, ξ₁)
            for ξ in range(ξ₁, 100, 10)
                @test abs(P(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ)
            end
        end
    end

    @testset "C_E" begin
        for (κ, λ) in params
            C = GinzburgLandauSelfSimilarSingular.C_E(κ, λ, ξ₁)
            for ξ in range(ξ₁, 100, 10)
                z = -im * κ / (1 - im * λ.ϵ) * ξ^2 / 2
                @test abs(GinzburgLandauSelfSimilarSingular.E(ξ, (λ, κ))) <=
                      C * exp(real(z)) * ξ^(-λ.d / 2 + 1 / λ.σ)
            end
        end
    end

    @testset "C_W" begin
        # This is not a proper test. It only checks the assumption for
        # W that is used for C_K.
        for (κ, λ) in params
            for ξ in range(Arb(1), 100, 10)
                d, ω, σ, ϵ = λ.d, λ.ω, λ.σ, λ.ϵ

                a = (1 / σ + im * ω / κ) / 2
                b = Arb(d) / 2
                z = -im * κ / (1 - im * ϵ) * ξ^2 / 2

                q1 = abs(GinzburgLandauSelfSimilarSingular.W(ξ, (λ, κ)))
                q2 =
                    κ / abs(1 - im * ϵ) *
                    exp(real(z)) *
                    exp(-imag(b - a) * π) *
                    abs(κ / 2(1 - im * ϵ))^-b *
                    ξ^(1 - d)

                inv_q1 = inv(q1)
                inv_q2 =
                    exp(-real(z)) *
                    abs(1 - im * ϵ) *
                    ξ^(d - 1) *
                    abs(κ / 2(1 - im * ϵ))^b *
                    exp(imag(b - a) * π) / κ

                @test Arblib.overlaps(q1, q2)
                @test Arblib.overlaps(inv_q1, inv_q2)
            end
        end
    end

    @testset "C_K" begin
        for (κ, λ) in params
            C = GinzburgLandauSelfSimilarSingular.C_K(κ, λ, ξ₁)
            for ξ in range(ξ₁, 100, 10)
                for η in range(ξ₁, 100, 10)
                    if η <= ξ
                        @test abs(GinzburgLandauSelfSimilarSingular.K(ξ, η, (λ, κ))) <=
                              C * ξ^(-1 / λ.σ) * η^(-1 + 1 / λ.σ)
                    else
                        @test abs(GinzburgLandauSelfSimilarSingular.K(ξ, η, (λ, κ))) <=
                              C * ξ^(-λ.d + 1 / λ.σ) * η^(-1 - 1 / λ.σ + λ.d)
                    end
                end
            end
        end
    end

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
