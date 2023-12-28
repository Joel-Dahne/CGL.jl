@testset "functions" begin
    ξ = Arb(30)
    κ = Arb(0.493223)
    λ = CGLParams{Arb}(1, 1.0, 2.3, 0.1, 0.2)
    _, _, _, c, c_dκ = CGL._abc_dκ(κ, λ)

    ξF64 = Float64(ξ)
    κF64 = Float64(κ)
    λF64 = CGLParams{Float64}(λ)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)
    fdm2 = central_fdm(5, 2)
    fdm3 = central_fdm(5, 3)

    @testset "P" begin
        @test Arblib.overlaps(P(ArbSeries((ξ, 1)), (λ, κ))[1], P_dξ(ξ, (λ, κ)))

        @test Arblib.overlaps(P(ξ, (λ, ArbSeries((κ, 1))))[1], P_dκ(ξ, (λ, κ)))

        @test CGL.P_asym_approx(ξ, (λ, κ)) ≈ P(ξ, (λ, κ)) rtol = 1e-20

        @test CGL.P_dξ_asym_approx(ξ, (λ, κ)) ≈ P_dξ(ξ, (λ, κ)) rtol = 1e-20

        @test P_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> P(ξ, (λF64, κF64)), ξF64) rtol = 1e-12

        @test P_dξ_dξ(ξ, (λ, κ)) ≈ fdm2(ξ -> P(ξ, (λF64, κF64)), ξF64) rtol = 1e-10

        @test P_dξ_dξ_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> P_dξ_dξ(ξ, (λF64, κF64)), ξF64) rtol = 1e-12
        @test P_dξ_dξ_dξ(ξ, (λ, κ)) ≈ fdm2(ξ -> P_dξ(ξ, (λF64, κF64)), ξF64) rtol = 1e-8
        @test P_dξ_dξ_dξ(ξ, (λ, κ)) ≈ fdm3(ξ -> P(ξ, (λF64, κF64)), ξF64) rtol = 1e-6

        @test P_dκ(ξ, (λ, κ)) ≈ fdm(κ -> P(ξF64, (λF64, κ)), κF64) rtol = 1e-12

        @test P_dξ_dκ(ξ, (λ, κ)) ≈ fdm(κ -> P_dξ(ξF64, (λF64, κ)), κF64) rtol = 1e-12

        @test P_dξ_dξ_dκ(ξ, (λ, κ)) ≈ fdm(κ -> P_dξ_dξ(ξF64, (λF64, κ)), κF64) rtol = 1e-12
    end

    @testset "E" begin
        @test Arblib.overlaps(E(ArbSeries((ξ, 1)), (λ, κ))[1], E_dξ(ξ, (λ, κ)))

        @test Arblib.overlaps(E(ξ, (λ, ArbSeries((κ, 1))))[1], E_dκ(ξ, (λ, κ)))

        @test E_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> E(ξ, (λF64, κF64)), ξF64) rtol = 1e-10

        @test E_dξ_dξ(ξ, (λ, κ)) ≈ fdm2(ξ -> E(ξ, (λF64, κF64)), ξF64) rtol = 1e-6

        @test E_dξ_dξ_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> E_dξ_dξ(ξ, (λF64, κF64)), ξF64) rtol = 1e-10
        @test E_dξ_dξ_dξ(ξ, (λ, κ)) ≈ fdm2(ξ -> E_dξ(ξ, (λF64, κF64)), ξF64) rtol = 1e-7
        @test E_dξ_dξ_dξ(ξ, (λ, κ)) ≈ fdm3(ξ -> E(ξ, (λF64, κF64)), ξF64) rtol = 1e-4

        @test E_dκ(ξ, (λ, κ)) ≈ fdm(κ -> E(ξF64, (λF64, κ)), κF64) rtol = 1e-10

        @test E_dξ_dκ(ξ, (λ, κ)) ≈ fdm(κ -> E_dξ(ξF64, (λF64, κ)), κF64) rtol = 1e-10
    end

    @testset "B_W" begin
        @test CGL.B_W_dκ(κ, λ) ≈ fdm(κ -> CGL.B_W(κ, λF64), κF64) rtol = 1e-12

        @test CGL.B_W_dκ(κ, λ) ≈ CGL.B_W_dκ(κF64, λF64) rtol = 1e-15
    end

    @testset "J_P" begin
        @test Arblib.overlaps(
            J_P(ξ, (λ, κ)),
            CGL.B_W(κ, λ) * P(ξ, (λ, κ)) * exp(-c * ξ^2) * ξ^(λ.d - 1),
        )

        @test J_P_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> J_P(ξ, (λF64, κF64)), ξF64) rtol = 1e-10

        @test J_P_dξ_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> J_P_dξ(ξ, (λF64, κF64)), ξF64) rtol = 1e-10
        @test J_P_dξ_dξ(ξ, (λ, κ)) ≈ fdm2(ξ -> J_P(ξ, (λF64, κF64)), ξF64) rtol = 1e-8

        @test J_P_dκ(ξ, (λ, κ)) ≈ fdm(κ -> J_P(ξF64, (λF64, κ)), κF64) rtol = 1e-10
    end

    @testset "J_E" begin
        @test Arblib.overlaps(
            J_E(ξ, (λ, κ)),
            CGL.B_W(κ, λ) * E(ξ, (λ, κ)) * exp(-c * ξ^2) * ξ^(λ.d - 1),
        )

        @test J_E_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> J_E(ξ, (λF64, κF64)), ξF64) rtol = 1e-10

        @test J_E_dξ_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> J_E_dξ(ξ, (λF64, κF64)), ξF64) rtol = 1e-10
        @test J_E_dξ_dξ(ξ, (λ, κ)) ≈ fdm2(ξ -> J_E(ξ, (λF64, κF64)), ξF64) rtol = 1e-8

        @test J_E_dκ(ξ, (λ, κ)) ≈ fdm(κ -> J_E(ξF64, (λF64, κ)), κF64) rtol = 1e-10
    end

    @testset "D" begin
        @test Arblib.overlaps(
            D(ξ, (λ, κ)),
            -c_dκ * CGL.B_W(κ, λ) * P(ξ, (λ, κ)) +
            CGL.B_W_dκ(κ, λ) * P(ξ, (λ, κ)) * ξ^-2 +
            CGL.B_W(κ, λ) * P_dκ(ξ, (λ, κ)) * ξ^-2,
        )

        @test D_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> D(ξ, (λF64, κF64)), ξF64) rtol = 1e-12

        @test D_dξ_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> D_dξ(ξ, (λF64, κF64)), ξF64) rtol = 1e-12
        @test D_dξ_dξ(ξ, (λ, κ)) ≈ fdm2(ξ -> D(ξ, (λF64, κF64)), ξF64) rtol = 1e-8
    end
end
