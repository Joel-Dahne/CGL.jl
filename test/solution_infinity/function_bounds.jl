@testset "solution_infinity_function_bounds" begin
    params = [
        (Arb(0.493223), Arb(0.0), CGLParams{Arb}(1, 1.0, 2.3, 0.0, 0.0)),
        (Arb(0.917383), Arb(0.0), CGLParams{Arb}(3, 1.0, 1.0, 0.0, 0.0)),
        (Arb(0.917383), Arb(0.01), CGLParams{Arb}(3, 1.0, 1.0, 0.01, 0.02)),
    ]

    ξ₁ = Arb(30)

    @testset "U $i" for (i, (κ, ϵ, λ)) in enumerate(params)
        (; d, σ) = λ
        a, b, c = CGL._abc(κ, ϵ, λ)
        z₁ = c * ξ₁^2

        CU = CGL.C_U(a, b, z₁)
        CU_dz = CGL.C_U_dz(a, b, z₁)
        CU_da = CGL.C_U_da(a, b, z₁)

        for z in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64] .* z₁
            @test abs(U(a, b, z)) <= CU * abs(z^(-a))
            @test abs(U(a, b, z)) >= 0.9CU * abs(z^(-a))

            @test abs(U_dz(a, b, z)) <= CU_dz * abs(z^(-a - 1))
            @test abs(U_dz(a, b, z)) >= 0.9CU_dz * abs(z^(-a - 1))

            @test abs(U_da(a, b, z)) <= CU_da * abs(log(z) * z^(-a))
            @test abs(U_da(a, b, z)) >= 0.9CU_da * abs(log(z) * z^(-a))
        end
    end

    @testset "FunctionBounds $i" for (i, (κ, ϵ, λ)) in enumerate(params)
        (; d, σ) = λ
        _, _, c = CGL._abc(κ, ϵ, λ)

        C = CGL.FunctionBounds(κ, ϵ, ξ₁, λ, include_dκ = true, include_dϵ = true)

        for ξ in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64] .* ξ₁
            ####
            ## P
            ####
            @test abs(P(ξ, κ, ϵ, λ)) <= C.P * ξ^(-1 / σ)
            @test abs(P(ξ, κ, ϵ, λ)) >= 0.9C.P * ξ^(-1 / σ)

            @test abs(P_dξ(ξ, κ, ϵ, λ)) <= C.P_dξ * ξ^(-1 / σ - 1)
            @test abs(P_dξ(ξ, κ, ϵ, λ)) >= 0.9C.P_dξ * ξ^(-1 / σ - 1)

            @test abs(P_dξ_dξ(ξ, κ, ϵ, λ)) <= C.P_dξ_dξ * ξ^(-1 / σ - 2)
            @test abs(P_dξ_dξ(ξ, κ, ϵ, λ)) >= 0.9C.P_dξ_dξ * ξ^(-1 / σ - 2)

            @test abs(P_dξ_dξ_dξ(ξ, κ, ϵ, λ)) <= C.P_dξ_dξ_dξ * ξ^(-1 / σ - 3)
            @test abs(P_dξ_dξ_dξ(ξ, κ, ϵ, λ)) >= 0.9C.P_dξ_dξ_dξ * ξ^(-1 / σ - 3)

            # IMPROVE: The three below ones don't give very tight bounds

            @test abs(P_dκ(ξ, κ, ϵ, λ)) <= C.P_dκ * log(ξ) * ξ^(-1 / σ)
            @test abs(P_dκ(ξ, κ, ϵ, λ)) >= 0.4C.P_dκ * log(ξ) * ξ^(-1 / σ)

            @test abs(P_dξ_dκ(ξ, κ, ϵ, λ)) <= C.P_dξ_dκ * log(ξ) * ξ^(-1 / σ - 1)
            @test abs(P_dξ_dκ(ξ, κ, ϵ, λ)) >= 0.25C.P_dξ_dκ * log(ξ) * ξ^(-1 / σ - 1)

            @test abs(P_dξ_dξ_dκ(ξ, κ, ϵ, λ)) <= C.P_dξ_dξ_dκ * log(ξ) * ξ^(-1 / σ - 2)
            @test abs(P_dξ_dξ_dκ(ξ, κ, ϵ, λ)) >= 0.1C.P_dξ_dξ_dκ * log(ξ) * ξ^(-1 / σ - 2)

            @test abs(P_dϵ(ξ, κ, ϵ, λ)) <= C.P_dϵ * ξ^(-1 / σ)
            @test abs(P_dϵ(ξ, κ, ϵ, λ)) >= 0.9C.P_dϵ * ξ^(-1 / σ)

            # IMPROVE: The two below ones don't give very tight bounds

            @test abs(P_dξ_dϵ(ξ, κ, ϵ, λ)) <= C.P_dξ_dϵ * ξ^(-1 / σ - 1)
            @test abs(P_dξ_dϵ(ξ, κ, ϵ, λ)) >= 0.25C.P_dξ_dϵ * ξ^(-1 / σ - 1)

            @test abs(P_dξ_dξ_dϵ(ξ, κ, ϵ, λ)) <= C.P_dξ_dξ_dϵ * ξ^(-1 / σ - 2)
            @test abs(P_dξ_dξ_dϵ(ξ, κ, ϵ, λ)) >= 0.05C.P_dξ_dξ_dϵ * ξ^(-1 / σ - 2)

            ####
            ## E
            ####
            @test abs(CGL.E(ξ, κ, ϵ, λ)) <= C.E * exp(real(c * ξ^2)) * ξ^(1 / σ - d)
            @test abs(CGL.E(ξ, κ, ϵ, λ)) >= 0.9C.E * exp(real(c * ξ^2)) * ξ^(1 / σ - d)

            @test abs(E_dξ(ξ, κ, ϵ, λ)) <= C.E_dξ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 1)
            @test abs(E_dξ(ξ, κ, ϵ, λ)) >=
                  0.9C.E_dξ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 1)

            @test abs(E_dξ_dξ(ξ, κ, ϵ, λ)) <=
                  C.E_dξ_dξ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 2)
            @test abs(E_dξ_dξ(ξ, κ, ϵ, λ)) >=
                  0.9C.E_dξ_dξ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 2)

            @test abs(E_dξ_dξ_dξ(ξ, κ, ϵ, λ)) <=
                  C.E_dξ_dξ_dξ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 3)
            @test abs(E_dξ_dξ_dξ(ξ, κ, ϵ, λ)) >=
                  0.9C.E_dξ_dξ_dξ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 3)

            @test abs(E_dκ(ξ, κ, ϵ, λ)) <= C.E_dκ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 2)
            @test abs(E_dκ(ξ, κ, ϵ, λ)) >=
                  0.9C.E_dκ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 2)

            @test abs(E_dξ_dκ(ξ, κ, ϵ, λ)) <=
                  C.E_dξ_dκ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 3)
            @test abs(E_dξ_dκ(ξ, κ, ϵ, λ)) >=
                  0.9C.E_dξ_dκ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 3)

            @test abs(E_dϵ(ξ, κ, ϵ, λ)) <= C.E_dϵ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 2)
            @test abs(E_dϵ(ξ, κ, ϵ, λ)) >=
                  0.9C.E_dϵ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 2)

            @test abs(E_dξ_dϵ(ξ, κ, ϵ, λ)) <=
                  C.E_dξ_dϵ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 3)
            @test abs(E_dξ_dϵ(ξ, κ, ϵ, λ)) >=
                  0.9C.E_dξ_dϵ * exp(real(c * ξ^2)) * ξ^(1 / σ - d + 3)

            ######
            ## J_P
            ######
            @test abs(J_P(ξ, κ, ϵ, λ)) <= C.J_P * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d - 1)
            @test abs(J_P(ξ, κ, ϵ, λ)) >=
                  0.9C.J_P * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d - 1)

            @test abs(J_P_dξ(ξ, κ, ϵ, λ)) <= C.J_P_dξ * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d)
            @test abs(J_P_dξ(ξ, κ, ϵ, λ)) >=
                  0.9C.J_P_dξ * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d)

            @test abs(J_P_dξ_dξ(ξ, κ, ϵ, λ)) <=
                  C.J_P_dξ_dξ * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d + 1)
            @test abs(J_P_dξ_dξ(ξ, κ, ϵ, λ)) >=
                  0.9C.J_P_dξ_dξ * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d + 1)

            @test abs(J_P_dκ(ξ, κ, ϵ, λ)) <=
                  C.J_P_dκ * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d + 1)
            @test abs(J_P_dκ(ξ, κ, ϵ, λ)) >=
                  0.9C.J_P_dκ * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d + 1)

            @test abs(J_P_dϵ(ξ, κ, ϵ, λ)) <=
                  C.J_P_dϵ * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d + 1)
            @test abs(J_P_dϵ(ξ, κ, ϵ, λ)) >=
                  0.9C.J_P_dϵ * exp(-real(c) * ξ^2) * ξ^(-1 / σ + d + 1)

            ######
            ## J_E
            ######
            @test abs(J_E(ξ, κ, ϵ, λ)) <= C.J_E * ξ^(1 / σ - 1)
            @test abs(J_E(ξ, κ, ϵ, λ)) >= 0.9C.J_E * ξ^(1 / σ - 1)

            @test abs(J_E_dξ(ξ, κ, ϵ, λ)) <= C.J_E_dξ * ξ^(1 / σ - 2)
            @test abs(J_E_dξ(ξ, κ, ϵ, λ)) >= 0.9C.J_E_dξ * ξ^(1 / σ - 2)

            @test abs(J_E_dξ_dξ(ξ, κ, ϵ, λ)) <= C.J_E_dξ_dξ * ξ^(1 / σ - 3)
            @test abs(J_E_dξ_dξ(ξ, κ, ϵ, λ)) >= 0.9C.J_E_dξ_dξ * ξ^(1 / σ - 3)

            @test abs(J_E_dκ(ξ, κ, ϵ, λ)) <= C.J_E_dκ * log(ξ) * ξ^(1 / σ - 1)
            # IMPROVE: This doesn't give a very tight bound
            @test abs(J_E_dκ(ξ, κ, ϵ, λ)) >= 0.3C.J_E_dκ * log(ξ) * ξ^(1 / σ - 1)

            @test abs(J_E_dϵ(ξ, κ, ϵ, λ)) <= C.J_E_dϵ * ξ^(1 / σ - 1)
            # IMPROVE: This doesn't give a very tight bound
            @test abs(J_E_dϵ(ξ, κ, ϵ, λ)) >= 0.2C.J_E_dϵ * ξ^(1 / σ - 1)

            ####
            ## D
            ####
            @test abs(D(ξ, κ, ϵ, λ)) <= C.D * ξ^(-1 / σ)
            @test abs(D(ξ, κ, ϵ, λ)) >= 0.9C.D * ξ^(-1 / σ)

            # IMPROVE: The two below ones don't give very tight bounds

            @test abs(D_dξ(ξ, κ, ϵ, λ)) <= C.D_dξ * ξ^(-1 / σ - 1)
            @test abs(D_dξ(ξ, κ, ϵ, λ)) >= 0.8C.D_dξ * ξ^(-1 / σ - 1)

            @test abs(D_dξ_dξ(ξ, κ, ϵ, λ)) <= C.D_dξ_dξ * ξ^(-1 / σ - 2)
            @test abs(D_dξ_dξ(ξ, κ, ϵ, λ)) >= 0.75C.D_dξ_dξ * ξ^(-1 / σ - 2)

            ####
            ## H
            ####
            @test abs(H(ξ, κ, ϵ, λ)) <= C.H * ξ^(-1 / σ)
            @test abs(H(ξ, κ, ϵ, λ)) >= 0.9C.H * ξ^(-1 / σ)

            # IMPROVE: The two below ones don't give very tight bounds

            @test abs(H_dξ(ξ, κ, ϵ, λ)) <= C.H_dξ * ξ^(-1 / σ - 1)
            @test abs(H_dξ(ξ, κ, ϵ, λ)) >= 0.9C.H_dξ * ξ^(-1 / σ - 1)

            @test abs(H_dξ_dξ(ξ, κ, ϵ, λ)) <= C.H_dξ_dξ * ξ^(-1 / σ - 2)
            @test abs(H_dξ_dξ(ξ, κ, ϵ, λ)) >= 0.9C.H_dξ_dξ * ξ^(-1 / σ - 2)
        end
    end
end
