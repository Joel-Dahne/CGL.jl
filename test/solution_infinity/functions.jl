@testset "functions" begin
    ξ = Arb(30)
    κ = Arb(0.493223)
    λ = gl_params(Arb, 1, 1.0, 2.3, 0.1, 0.2)
    _, _, c = CGL._abc(κ, λ)

    ξF64 = Float64(ξ)
    κF64 = Float64(κ)
    λF64 = gl_params(Float64, λ)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)
    fdm2 = central_fdm(5, 2)

    @testset "P" begin
        @test Arblib.overlaps(P(ArbSeries((ξ, 1)), (λ, κ))[1], P_dξ(ξ, (λ, κ)))

        @test Arblib.overlaps(P(ξ, (λ, ArbSeries((κ, 1))))[1], P_dκ(ξ, (λ, κ)))

        @test CGL.P_asym_approx(ξ, (λ, κ)) ≈ P(ξ, (λ, κ)) rtol = 1e-20

        @test CGL.P_dξ_asym_approx(ξ, (λ, κ)) ≈ P_dξ(ξ, (λ, κ)) rtol = 1e-20

        @test P_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> P(ξ, (λF64, κF64)), ξF64) rtol = 1e-12

        @test P_dξ_dξ(ξ, (λ, κ)) ≈ fdm2(ξ -> P(ξ, (λF64, κF64)), ξF64) rtol = 1e-10

        @test P_dκ(ξ, (λ, κ)) ≈ fdm(κ -> P(ξF64, (λF64, κ)), κF64) rtol = 1e-12

        @test P_dξ_dκ(ξ, (λ, κ)) ≈ fdm(κ -> P_dξ(ξF64, (λF64, κ)), κF64) rtol = 1e-12
    end

    @testset "E" begin
        @test Arblib.overlaps(E(ArbSeries((ξ, 1)), (λ, κ))[1], E_dξ(ξ, (λ, κ)))

        @test Arblib.overlaps(E(ξ, (λ, ArbSeries((κ, 1))))[1], E_dκ(ξ, (λ, κ)))

        @test E_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> E(ξ, (λF64, κF64)), ξF64) rtol = 1e-10

        @test E_dξ_dξ(ξ, (λ, κ)) ≈ fdm2(ξ -> E(ξ, (λF64, κF64)), ξF64) rtol = 1e-6

        @test E_dκ(ξ, (λ, κ)) ≈ fdm(κ -> E(ξF64, (λF64, κ)), κF64) rtol = 1e-10

        @test E_dξ_dκ(ξ, (λ, κ)) ≈ fdm(κ -> E_dξ(ξF64, (λF64, κ)), κF64) rtol = 1e-10
    end

    @testset "B_W" begin
        @test CGL.B_W_dκ(κ, λ) ≈ fdm(κ -> CGL.B_W(κ, λF64), κF64) rtol = 1e-12
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

        @test J_P_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> J_P(ξ, (λF64, κF64)), ξF64) rtol = 1e-10

        @test J_P_dκ(ξ, (λ, κ)) ≈ fdm(κ -> J_P(ξF64, (λF64, κ)), κF64) rtol = 1e-10
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

        @test J_E_dξ(ξ, (λ, κ)) ≈ fdm(ξ -> J_E(ξ, (λF64, κF64)), ξF64) rtol = 1e-10

        @test J_E_dκ(ξ, (λ, κ)) ≈ fdm(κ -> J_E(ξF64, (λF64, κ)), κF64) rtol = 1e-10
    end
end
