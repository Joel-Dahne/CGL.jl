@testset "solution_infinity_function_bounds" begin
    params = [
        (Arb(0.493223), gl_params(Arb, 1, 1.0, 2.3, 0.0, 0.0)),
        (Arb(0.917383), gl_params(Arb, 3, 1.0, 1.0, 0.0, 0.0)),
        (Arb(0.917383), gl_params(Arb, 3, 1.0, 1.0, 0.01, 0.02)),
    ]

    ξ₁ = Arb(10)

    @testset "C_hypgeom_u" begin
        for (κ, λ) in params
            a, b, c = CGL._abc(κ, λ)
            z₁ = c * ξ₁^2
            C = CGL.C_hypgeom_u(a, b, z₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                z = k * z₁
                @test abs(hypgeom_u(a, b, z)) <= C * abs(z^(-a))
            end
        end
    end

    @testset "C_hypgeom_u_dz" begin
        for (κ, λ) in params
            a, b, c = CGL._abc(κ, λ)
            z₁ = c * ξ₁^2

            C = CGL.C_hypgeom_u_dz(a, b, z₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                z = k * z₁
                @test abs(hypgeom_u_dz(a, b, z)) <= C * abs(z^(-a))
            end
        end
    end

    @testset "C_P" begin
        for (κ, λ) in params
            C = CGL.C_P(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(P(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ)
            end
        end
    end

    @testset "C_E" begin
        for (κ, λ) in params
            _, _, c = CGL._abc(κ, λ)
            C = CGL.C_E(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(CGL.E(ξ, (λ, κ))) <= C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d)
            end
        end
    end

    @testset "C_P_dξ" begin
        for (κ, λ) in params
            C = CGL.C_P_dξ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(P_dξ(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ - 1)
            end
        end
    end

    @testset "C_E_dξ" begin
        for (κ, λ) in params
            _, _, c = CGL._abc(κ, λ)
            C = CGL.C_E_dξ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(E_dξ(ξ, (λ, κ))) <= C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 1)
            end
        end
    end

    @testset "C_P_dξ_dξ" begin
        for (κ, λ) in params
            C = CGL.C_P_dξ_dξ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(P_dξ_dξ(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ - 2)
            end
        end
    end

    @testset "C_P_dξ_dξ" begin
        for (κ, λ) in params
            C = CGL.C_P_dξ_dξ_dξ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(P_dξ_dξ_dξ(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ - 3)
            end
        end
    end

    @testset "C_E_dξ_dξ" begin
        for (κ, λ) in params
            _, _, c = CGL._abc(κ, λ)
            C = CGL.C_E_dξ_dξ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(E_dξ_dξ(ξ, (λ, κ))) <=
                      C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 2)
            end
        end
    end

    @testset "C_P_dκ" begin
        for (κ, λ) in params
            C = CGL.C_P_dκ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(P_dκ(ξ, (λ, κ))) <= C * log(ξ) * ξ^(-1 / λ.σ)
            end
        end
    end

    @testset "C_E_dκ" begin
        for (κ, λ) in params
            _, _, c = CGL._abc(κ, λ)
            C = CGL.C_E_dκ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(E_dκ(ξ, (λ, κ))) <= C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 2)
            end
        end
    end

    @testset "C_P_dξ_dκ" begin
        # The implementation of C_P_dξ_dκ gives terrible bounds when ξ
        # is to small. We therefore take a large ξ in this case to be
        # able to check anything. The part giving bad bounds is
        # _hypgeom_u_da_finite_difference
        let ξ₁ = 2ξ₁
            for (κ, λ) in params
                C = CGL.C_P_dξ_dκ(κ, λ, ξ₁)
                for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                    ξ = k * ξ₁
                    @test abs(P_dξ_dκ(ξ, (λ, κ))) <= C * log(ξ) * ξ^(-1 / λ.σ - 1)
                end
            end
        end
    end

    @testset "C_E_dξ_dκ" begin
        # The implementation of C_E_dξ_dκ gives terrible bounds when ξ
        # is to small. We therefore take a large ξ in this case to be
        # able to check anything. The part giving bad bounds is
        # _hypgeom_u_da_finite_difference
        let ξ₁ = 2ξ₁
            for (κ, λ) in params
                _, _, c = CGL._abc(κ, λ)
                C = CGL.C_E_dξ_dκ(κ, λ, ξ₁)
                for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                    ξ = k * ξ₁
                    @test abs(E_dξ_dκ(ξ, (λ, κ))) <=
                          C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 3)
                end
            end
        end
    end

    @testset "C_W" begin
        # This is not a proper test. It only checks the assumption for
        # W that is used for C_K.
        for (κ, λ) in params
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                d, ω, σ, ϵ = λ.d, λ.ω, λ.σ, λ.ϵ

                a = (1 / σ + im * ω / κ) / 2
                b = Arb(d) / 2
                z = -im * κ / (1 - im * ϵ) * ξ^2 / 2

                sgn = sign(imag(z))

                q1 = abs(CGL.W(ξ, (λ, κ)))
                q2 =
                    κ / abs(1 - im * ϵ) *
                    exp(real(z)) *
                    exp(-sgn * imag(b - a) * π) *
                    abs(κ / 2(1 - im * ϵ))^-b *
                    ξ^(1 - d)

                inv_q1 = inv(q1)
                inv_q2 =
                    exp(-real(z)) *
                    abs(1 - im * ϵ) *
                    ξ^(d - 1) *
                    abs(κ / 2(1 - im * ϵ))^b *
                    exp(sgn * imag(b - a) * π) / κ

                @test Arblib.overlaps(q1, q2)
                @test Arblib.overlaps(inv_q1, inv_q2)
            end
        end
    end

    @testset "C_K" begin
        for (κ, λ) in params
            C = CGL.C_K(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                for l in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                    η = l * ξ₁
                    if η <= ξ
                        @test abs(CGL.K(ξ, η, (λ, κ))) <=
                              C * ξ^(-1 / λ.σ) * η^(-1 + 1 / λ.σ)
                    else
                        @test abs(CGL.K(ξ, η, (λ, κ))) <=
                              C * ξ^(-λ.d + 1 / λ.σ) * η^(-1 - 1 / λ.σ + λ.d)
                    end
                end
            end
        end
    end

    @testset "C_J_P" begin
        for (κ, λ) in params
            _, _, c = CGL._abc(κ, λ)
            C = CGL.C_J_P(κ, ξ₁, λ)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(J_P(ξ, (λ, κ))) <=
                      C * exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d - 1)
            end
        end
    end

    @testset "C_J_E" begin
        for (κ, λ) in params
            C = CGL.C_J_E(κ, ξ₁, λ)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(J_E(ξ, (λ, κ))) <= C * ξ^(1 / λ.σ - 1)
            end
        end
    end

    @testset "C_J_P_dξ" begin
        for (κ, λ) in params
            _, _, c = CGL._abc(κ, λ)
            C = CGL.C_J_P_dξ(κ, ξ₁, λ)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(J_P_dξ(ξ, (λ, κ))) <= C * exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d)
            end
        end
    end

    @testset "C_J_E_dξ" begin
        for (κ, λ) in params
            C = CGL.C_J_E_dξ(κ, ξ₁, λ)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(J_E_dξ(ξ, (λ, κ))) <= C * ξ^(1 / λ.σ - 2)
            end
        end
    end

    @testset "C_J_P_dκ" begin
        for (κ, λ) in params
            _, _, c = CGL._abc(κ, λ)
            C = CGL.C_J_P_dκ(κ, ξ₁, λ)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(J_P_dκ(ξ, (λ, κ))) <=
                      C * exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d + 1)
            end
        end
    end

    @testset "C_J_E_dκ" begin
        # The implementation of C_P_dξ_dκ gives terrible bounds when ξ
        # is to small. We therefore take a large ξ in this case to be
        # able to check anything. The part giving bad bounds is
        # _hypgeom_u_da_finite_difference
        let ξ₁ = 2ξ₁
            for (κ, λ) in params
                C = CGL.C_J_E_dκ(κ, ξ₁, λ)
                for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                    ξ = k * ξ₁
                    @test abs(J_E_dκ(ξ, (λ, κ))) <= C * log(ξ) * ξ^(1 / λ.σ - 1)
                end
            end
        end
    end

    @testset "D" begin
        for (κ, λ) in params
            C = CGL.C_D(κ, ξ₁, λ)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(D(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ)
            end
        end
    end

    @testset "D_dξ" begin
        # The implementation of C_D_dξ gives terrible bounds when ξ is
        # to small. We therefore take a large ξ in this case to be
        # able to check anything. The part giving bad bounds is
        # _hypgeom_u_da_finite_difference
        let ξ₁ = 2ξ₁
            for (κ, λ) in params
                C = CGL.C_D_dξ(κ, ξ₁, λ)
                for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                    ξ = k * ξ₁
                    @test abs(D_dξ(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ - 1)
                end
            end
        end
    end

    @testset "D_dξ_dξ" begin
        # The implementation of C_D_dξ_dξ gives terrible bounds when ξ
        # is to small. We therefore take a large ξ in this case to be
        # able to check anything. The part giving bad bounds is
        # _hypgeom_u_da_finite_difference
        let ξ₁ = 2ξ₁
            for (κ, λ) in params
                C = CGL.C_D_dξ_dξ(κ, ξ₁, λ)
                for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                    ξ = k * ξ₁
                    @test abs(D_dξ_dξ(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ - 2)
                end
            end
        end
    end
end
