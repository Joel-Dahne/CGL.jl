@testset "functions" begin
    ξ = Arb(30)
    κ = Arb(0.493223)
    ϵ = Arb(0.1)
    Λ = CGLParams{Arb}(1, 1.0, 2.3, 0.2)
    _, _, _, c, c_dκ = CGL._abc_dκ(κ, ϵ, Λ)
    _, _, c, c_dϵ = CGL._abc_dϵ(κ, ϵ, Λ)

    ξF64 = Float64(ξ)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    ΛF64 = CGLParams{Float64}(Λ)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)
    fdm2 = central_fdm(5, 2)
    fdm3 = central_fdm(5, 3)

    @testset "P" begin
        @test Arblib.overlaps(P(ArbSeries((ξ, 1)), κ, ϵ, Λ)[1], P_dξ(ξ, κ, ϵ, Λ))

        @test Arblib.overlaps(P(ξ, ArbSeries((κ, 1)), ϵ, Λ)[1], P_dκ(ξ, κ, ϵ, Λ))

        @test P_dξ(ξ, κ, ϵ, Λ) ≈ fdm(ξ -> P(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-12

        @test P_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm2(ξ -> P(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-10

        @test P_dξ_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm(ξ -> P_dξ_dξ(ξ, κF64, ϵF64, ΛF64), ξF64) rtol =
            1e-12
        @test P_dξ_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm2(ξ -> P_dξ(ξ, κF64, ϵF64, ΛF64), ξF64) rtol =
            1e-8
        @test P_dξ_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm3(ξ -> P(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-6

        @test P_dκ(ξ, κ, ϵ, Λ) ≈ fdm(κ -> P(ξF64, κ, ϵF64, ΛF64), κF64) rtol = 1e-12

        @test P_dξ_dκ(ξ, κ, ϵ, Λ) ≈ fdm(κ -> P_dξ(ξF64, κ, ϵF64, ΛF64), κF64) rtol = 1e-12

        @test P_dξ_dξ_dκ(ξ, κ, ϵ, Λ) ≈ fdm(κ -> P_dξ_dξ(ξF64, κ, ϵF64, ΛF64), κF64) rtol =
            1e-12

        @test P_dϵ(ξ, κ, ϵ, Λ) ≈ fdm(ϵ -> P(ξF64, κF64, ϵ, ΛF64), ϵF64) rtol = 1e-12

        @test P_dξ_dϵ(ξ, κ, ϵ, Λ) ≈ fdm(ϵ -> P_dξ(ξF64, κF64, ϵ, ΛF64), ϵF64) rtol = 1e-12

        @test P_dξ_dξ_dϵ(ξ, κ, ϵ, Λ) ≈ fdm(ϵ -> P_dξ_dξ(ξF64, κF64, ϵ, ΛF64), ϵF64) rtol =
            1e-12
    end

    @testset "E" begin
        @test Arblib.overlaps(E(ArbSeries((ξ, 1)), κ, ϵ, Λ)[1], E_dξ(ξ, κ, ϵ, Λ))

        @test Arblib.overlaps(E(ξ, ArbSeries((κ, 1)), ϵ, Λ)[1], E_dκ(ξ, κ, ϵ, Λ))

        @test E_dξ(ξ, κ, ϵ, Λ) ≈ fdm(ξ -> E(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-10

        @test E_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm2(ξ -> E(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-6

        @test E_dξ_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm(ξ -> E_dξ_dξ(ξ, κF64, ϵF64, ΛF64), ξF64) rtol =
            1e-10
        @test E_dξ_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm2(ξ -> E_dξ(ξ, κF64, ϵF64, ΛF64), ξF64) rtol =
            1e-7
        @test E_dξ_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm3(ξ -> E(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-4

        @test E_dκ(ξ, κ, ϵ, Λ) ≈ fdm(κ -> E(ξF64, κ, ϵF64, ΛF64), κF64) rtol = 1e-10

        @test E_dξ_dκ(ξ, κ, ϵ, Λ) ≈ fdm(κ -> E_dξ(ξF64, κ, ϵF64, ΛF64), κF64) rtol = 1e-10

        @test E_dϵ(ξ, κ, ϵ, Λ) ≈ fdm(ϵ -> E(ξF64, κF64, ϵ, ΛF64), ϵF64) rtol = 1e-10

        @test E_dξ_dϵ(ξ, κ, ϵ, Λ) ≈ fdm(ϵ -> E_dξ(ξF64, κF64, ϵ, ΛF64), ϵF64) rtol = 1e-10
    end

    @testset "B_W" begin
        @test CGL.B_W_dκ(κ, ϵ, Λ) ≈ fdm(κ -> CGL.B_W(κ, ϵF64, ΛF64), κF64) rtol = 1e-12

        @test CGL.B_W_dκ(κ, ϵ, Λ) ≈ CGL.B_W_dκ(κF64, ϵF64, ΛF64) rtol = 1e-15

        @test CGL.B_W_dϵ(κ, ϵ, Λ) ≈ fdm(ϵ -> CGL.B_W(κF64, ϵ, ΛF64), ϵF64) rtol = 1e-12
    end

    @testset "J_P" begin
        @test Arblib.overlaps(
            J_P(ξ, κ, ϵ, Λ),
            CGL.B_W(κ, ϵ, Λ) * P(ξ, κ, ϵ, Λ) * exp(-c * ξ^2) * ξ^(Λ.d - 1),
        )

        @test J_P_dκ(ξ, κ, ϵ, Λ) ≈ fdm(κ -> J_P(ξF64, κ, ϵF64, ΛF64), κF64) rtol = 1e-10

        @test J_P_dϵ(ξ, κ, ϵ, Λ) ≈ fdm(ϵ -> J_P(ξF64, κF64, ϵ, ΛF64), ϵF64) rtol = 1e-10
    end

    @testset "J_E" begin
        @test Arblib.overlaps(
            J_E(ξ, κ, ϵ, Λ),
            CGL.B_W(κ, ϵ, Λ) * E(ξ, κ, ϵ, Λ) * exp(-c * ξ^2) * ξ^(Λ.d - 1),
        )

        @test J_E_dκ(ξ, κ, ϵ, Λ) ≈ fdm(κ -> J_E(ξF64, κ, ϵF64, ΛF64), κF64) rtol = 1e-10

        @test J_E_dϵ(ξ, κ, ϵ, Λ) ≈ fdm(ϵ -> J_E(ξF64, κF64, ϵ, ΛF64), ϵF64) rtol = 1e-10
    end

    @testset "D" begin
        @test Arblib.overlaps(
            D(ξ, κ, ϵ, Λ),
            -c_dκ * CGL.B_W(κ, ϵ, Λ) * P(ξ, κ, ϵ, Λ) +
            CGL.B_W_dκ(κ, ϵ, Λ) * P(ξ, κ, ϵ, Λ) * ξ^-2 +
            CGL.B_W(κ, ϵ, Λ) * P_dκ(ξ, κ, ϵ, Λ) * ξ^-2,
        )

        @test D_dξ(ξ, κ, ϵ, Λ) ≈ fdm(ξ -> D(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-12

        @test D_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm(ξ -> D_dξ(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-12
        @test D_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm2(ξ -> D(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-8
    end

    @testset "H" begin
        @test Arblib.overlaps(
            H(ξ, κ, ϵ, Λ),
            -c_dϵ * CGL.B_W(κ, ϵ, Λ) * P(ξ, κ, ϵ, Λ) +
            CGL.B_W_dϵ(κ, ϵ, Λ) * P(ξ, κ, ϵ, Λ) * ξ^-2 +
            CGL.B_W(κ, ϵ, Λ) * P_dϵ(ξ, κ, ϵ, Λ) * ξ^-2,
        )

        @test H_dξ(ξ, κ, ϵ, Λ) ≈ fdm(ξ -> H(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-12

        @test H_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm(ξ -> H_dξ(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-12
        @test H_dξ_dξ(ξ, κ, ϵ, Λ) ≈ fdm2(ξ -> H(ξ, κF64, ϵF64, ΛF64), ξF64) rtol = 1e-8
    end
end
