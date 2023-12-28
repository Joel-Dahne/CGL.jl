@testset "solution_infinity_function_bounds" begin
    params = [
        (Arb(0.493223), CGLParams{Arb}(1, 1.0, 2.3, 0.0, 0.0)),
        (Arb(0.917383), CGLParams{Arb}(3, 1.0, 1.0, 0.0, 0.0)),
        (Arb(0.917383), CGLParams{Arb}(3, 1.0, 1.0, 0.01, 0.02)),
    ]

    ξ₁ = Arb(30)

    @testset "C_hypgeom_u" begin
        for (κ, λ) in params
            a, b, c = CGL._abc(κ, λ)
            z₁ = c * ξ₁^2
            C = CGL.C_hypgeom_u(a, b, z₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                z = k * z₁
                @test abs(hypgeom_u(a, b, z)) <= C * abs(z^(-a))
                @test abs(hypgeom_u(a, b, z)) >= 0.9C * abs(z^(-a))
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
                @test abs(hypgeom_u_dz(a, b, z)) <= C * abs(z^(-a - 1))
                @test abs(hypgeom_u_dz(a, b, z)) >= 0.9C * abs(z^(-a - 1))
            end
        end
    end

    @testset "C_hypgeom_u_da" begin
        for (κ, λ) in params
            a, b, c = CGL._abc(κ, λ)
            z₁ = c * ξ₁^2
            C = CGL.C_hypgeom_u_da(a, b, z₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                z = k * z₁
                @test abs(hypgeom_u_da(a, b, z)) <= C * abs(log(z) * z^(-a))
                @test abs(hypgeom_u_da(a, b, z)) >= 0.9C * abs(log(z) * z^(-a))
            end
        end
    end

    @testset "C_P" begin
        for (κ, λ) in params
            C = CGL.C_P(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(P(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ)
                @test abs(P(ξ, (λ, κ))) >= 0.9C * ξ^(-1 / λ.σ)
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
                @test abs(CGL.E(ξ, (λ, κ))) >= 0.9C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d)
            end
        end
    end

    @testset "C_P_dξ" begin
        for (κ, λ) in params
            C = CGL.C_P_dξ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(P_dξ(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ - 1)
                @test abs(P_dξ(ξ, (λ, κ))) >= 0.9C * ξ^(-1 / λ.σ - 1)
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
                @test abs(E_dξ(ξ, (λ, κ))) >=
                      0.9C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 1)
            end
        end
    end

    @testset "C_P_dξ_dξ" begin
        for (κ, λ) in params
            C = CGL.C_P_dξ_dξ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(P_dξ_dξ(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ - 2)
                @test abs(P_dξ_dξ(ξ, (λ, κ))) >= 0.9C * ξ^(-1 / λ.σ - 2)
            end
        end
    end

    @testset "C_P_dξ_dξ_dξ" begin
        for (κ, λ) in params
            C = CGL.C_P_dξ_dξ_dξ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(P_dξ_dξ_dξ(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ - 3)
                @test abs(P_dξ_dξ_dξ(ξ, (λ, κ))) >= 0.9C * ξ^(-1 / λ.σ - 3)
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
                @test abs(E_dξ_dξ(ξ, (λ, κ))) >=
                      0.9C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 2)
            end
        end
    end

    @testset "C_E_dξ_dξ_dξ" begin
        for (κ, λ) in params
            _, _, c = CGL._abc(κ, λ)
            C = CGL.C_E_dξ_dξ_dξ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(E_dξ_dξ_dξ(ξ, (λ, κ))) <=
                      C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 3)
                @test abs(E_dξ_dξ_dξ(ξ, (λ, κ))) >=
                      0.9C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 3)
            end
        end
    end

    @testset "C_P_dκ" begin
        for (κ, λ) in params
            C = CGL.C_P_dκ(κ, λ, ξ₁)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(P_dκ(ξ, (λ, κ))) <= C * log(ξ) * ξ^(-1 / λ.σ)
                # IMPROVE: This doesn't give a very tight upper bound
                @test abs(P_dκ(ξ, (λ, κ))) >= 0.4C * log(ξ) * ξ^(-1 / λ.σ)
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
                @test abs(E_dκ(ξ, (λ, κ))) >=
                      0.9C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 2)
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
                    # IMPROVE: This doesn't give a very tight upper bound
                    @test abs(P_dξ_dκ(ξ, (λ, κ))) >= 0.3C * log(ξ) * ξ^(-1 / λ.σ - 1)
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
                    @test abs(E_dξ_dκ(ξ, (λ, κ))) >=
                          0.9C * exp(real(c * ξ^2)) * ξ^(1 / λ.σ - λ.d + 3)
                end
            end
        end
    end

    @testset "C_P_dξ_dξ_dκ" begin
        # The implementation of C_P_dξ_dκ gives terrible bounds when ξ
        # is to small. We therefore take a large ξ in this case to be
        # able to check anything. The part giving bad bounds is
        # _hypgeom_u_da_finite_difference
        let ξ₁ = 2ξ₁
            for (κ, λ) in params
                C = CGL.C_P_dξ_dξ_dκ(κ, λ, ξ₁)
                for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                    ξ = k * ξ₁
                    @test abs(P_dξ_dξ_dκ(ξ, (λ, κ))) <= C * log(ξ) * ξ^(-1 / λ.σ - 2)
                    # IMPROVE: This doesn't give a very tight upper bound
                    @test abs(P_dξ_dξ_dκ(ξ, (λ, κ))) >= 0.1C * log(ξ) * ξ^(-1 / λ.σ - 2)
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
                @test abs(J_P(ξ, (λ, κ))) >=
                      0.9C * exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d - 1)
            end
        end
    end

    @testset "C_J_E" begin
        for (κ, λ) in params
            C = CGL.C_J_E(κ, ξ₁, λ)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(J_E(ξ, (λ, κ))) <= C * ξ^(1 / λ.σ - 1)
                @test abs(J_E(ξ, (λ, κ))) >= 0.9C * ξ^(1 / λ.σ - 1)
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
                @test abs(J_P_dξ(ξ, (λ, κ))) >=
                      0.9C * exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d)
            end
        end
    end

    @testset "C_J_E_dξ" begin
        for (κ, λ) in params
            C = CGL.C_J_E_dξ(κ, ξ₁, λ)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(J_E_dξ(ξ, (λ, κ))) <= C * ξ^(1 / λ.σ - 2)
                @test abs(J_E_dξ(ξ, (λ, κ))) >= 0.9C * ξ^(1 / λ.σ - 2)
            end
        end
    end

    @testset "C_J_P_dξ_dξ" begin
        for (κ, λ) in params
            _, _, c = CGL._abc(κ, λ)
            C = CGL.C_J_P_dξ_dξ(κ, ξ₁, λ)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(J_P_dξ_dξ(ξ, (λ, κ))) <=
                      C * exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d + 1)
                @test abs(J_P_dξ_dξ(ξ, (λ, κ))) >=
                      0.9C * exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d + 1)
            end
        end
    end

    @testset "C_J_E_dξ_dξ" begin
        for (κ, λ) in params
            C = CGL.C_J_E_dξ_dξ(κ, ξ₁, λ)
            for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                ξ = k * ξ₁
                @test abs(J_E_dξ_dξ(ξ, (λ, κ))) <= C * ξ^(1 / λ.σ - 3)
                @test abs(J_E_dξ_dξ(ξ, (λ, κ))) >= 0.9C * ξ^(1 / λ.σ - 3)
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
                @test abs(J_P_dκ(ξ, (λ, κ))) >=
                      0.9C * exp(-real(c) * ξ^2) * ξ^(-1 / λ.σ + λ.d + 1)
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
                    # IMPROVE: This doesn't give a very tight upper bound
                    @test abs(J_E_dκ(ξ, (λ, κ))) >= 0.3C * log(ξ) * ξ^(1 / λ.σ - 1)
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
                @test abs(D(ξ, (λ, κ))) >= 0.9C * ξ^(-1 / λ.σ)
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
                    @test abs(D_dξ(ξ, (λ, κ))) >= 0.9C * ξ^(-1 / λ.σ - 1)
                end
            end
        end
    end

    @testset "D_dξ_dξ" begin
        # The implementation of C_D_dξ_dξ gives terrible bounds when ξ
        # is to small. We therefore take a large ξ in this case to be
        # able to check anything. The part giving bad bounds is
        # _hypgeom_u_da_finite_difference
        let ξ₁ = 3ξ₁
            for (κ, λ) in params
                C = CGL.C_D_dξ_dξ(κ, ξ₁, λ)
                for k in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64]
                    ξ = k * ξ₁
                    @test abs(D_dξ_dξ(ξ, (λ, κ))) <= C * ξ^(-1 / λ.σ - 2)
                    @test abs(D_dξ_dξ(ξ, (λ, κ))) >= 0.9C * ξ^(-1 / λ.σ - 2)
                end
            end
        end
    end
end
