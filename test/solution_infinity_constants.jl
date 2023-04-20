@testset "solution_infinity_constants" begin
    paramss = [
        gl_params(1, Arb(1.0), 2.3, 0.0, 0.0),
        gl_params(3, Arb(1.41727), 1.0, 0.0, 0.0),
        gl_params(3, Arb(1.41727), 1.0, 0.01, 0.02),
    ]

    κs = Arb[0.49323, 0.45535, 0.45535]

    @testset "C_P" begin
        for (κ, p) in zip(κs, paramss)
            C = GinzburgLandauSelfSimilarSingular.C_P(κ, p)
            for ξ in range(Arb(1), 100, 10)
                @test abs(P(ξ, (p, κ))) <= C * ξ^(-1 / p.σ)
            end
        end
    end

    @testset "C_P" begin
        for (κ, p) in zip(κs, paramss)
            C = GinzburgLandauSelfSimilarSingular.C_E(κ, p)
            for ξ in range(Arb(1), 100, 10)
                z = -im * κ / (1 - im * p.ϵ) * ξ^2 / 2
                @test abs(GinzburgLandauSelfSimilarSingular.E(ξ, (p, κ))) <=
                      C * exp(real(z)) * ξ^(-p.d / 2 + 1 / p.σ)
            end
        end
    end

    @testset "C_W" begin
        # This is not a proper test. It only checks the assumption for
        # W that is used for C_K.
        for (κ, p) in zip(κs, paramss)
            for ξ in range(Arb(1), 100, 10)
                d, ω, σ, ϵ = p.d, p.ω, p.σ, p.ϵ

                a = (1 / σ + im * ω / κ) / 2
                b = Arb(d) / 2
                z = -im * κ / (1 - im * ϵ) * ξ^2 / 2

                q1 = abs(GinzburgLandauSelfSimilarSingular.W(ξ, (p, κ)))
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
        for (κ, p) in zip(κs, paramss)
            C = GinzburgLandauSelfSimilarSingular.C_K(κ, p)
            for ξ in range(Arb(1), 100, 10)
                for η in range(Arb(1), 100, 10)
                    if η <= ξ
                        @test abs(GinzburgLandauSelfSimilarSingular.K(ξ, η, (p, κ))) <=
                              C * ξ^(-1 / p.σ) * η^(-1 + 1 / p.σ)
                    else
                        @test abs(GinzburgLandauSelfSimilarSingular.K(ξ, η, (p, κ))) <=
                              C * ξ^(-p.d + 1 / p.σ) * η^(-1 - 1 / p.σ + p.d)
                    end
                end
            end
        end
    end

    @testset "C_T1" begin
        # We don't have an obvious way to check this inequality. For
        # now we simply check that C₁ is finite and C₂ it is
        # increasing in v.
        for (κ, p) in zip(κs, paramss)
            C1₁, C1₂ = GinzburgLandauSelfSimilarSingular.C_T1(Arb(0.1), κ, p)
            C2₁, C2₂ = GinzburgLandauSelfSimilarSingular.C_T1(Arb(0.2), κ, p)

            @test isfinite(C1₁)
            @test isfinite(C2₁)
            @test C1₂ < C2₂
        end
    end
end
