@testset "U_expansion" begin
    params = [
        (Arb(0.493223), Arb(0.0), CGLParams{Arb}(1, 1.0, 2.3, 0.0)),
        (Arb(0.917383), Arb(0.0), CGLParams{Arb}(3, 1.0, 1.0, 0.0)),
        (Arb(0.917383), Arb(0.01), CGLParams{Arb}(3, 1.0, 1.0, 0.02)),
    ]

    ξ₁ = Arb(30)

    @testset "U $i" for (i, (κ, ϵ, λ)) in enumerate(params)
        (; d, σ) = λ
        a, b, c = CGL._abc(κ, ϵ, λ)
        z₁ = c * ξ₁^2

        CU = CGL.UBounds(a, b, c, ξ₁, include_da = true)

        for z in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64] .* z₁
            # U

            @test abs(CGL.U(a, b, z)) <= CU.U_a_b * abs(z^(-a))
            @test abs(CGL.U(a, b, z)) >= 0.9CU.U_a_b * abs(z^(-a))

            @test abs(CGL.U(a + 1, b + 1, z)) <= CU.U_ap1_bp1 * abs(z^(-(a + 1)))
            @test abs(CGL.U(a + 1, b + 1, z)) >= 0.9CU.U_ap1_bp1 * abs(z^(-(a + 1)))

            @test abs(CGL.U(a + 2, b + 2, z)) <= CU.U_ap2_bp2 * abs(z^(-(a + 2)))
            @test abs(CGL.U(a + 2, b + 2, z)) >= 0.9CU.U_ap2_bp2 * abs(z^(-(a + 2)))

            @test abs(CGL.U(b - a, b, z)) <= CU.U_bma_b * abs(z^(-(b - a)))
            @test abs(CGL.U(b - a, b, z)) >= 0.9CU.U_bma_b * abs(z^(-(b - a)))

            @test abs(CGL.U(b - a + 1, b + 1, z)) <= CU.U_bmap1_bp1 * abs(z^(-(b - a + 1)))
            @test abs(CGL.U(b - a + 1, b + 1, z)) >=
                  0.9CU.U_bmap1_bp1 * abs(z^(-(b - a + 1)))

            @test abs(CGL.U(b - a + 2, b + 2, z)) <= CU.U_bmap2_bp2 * abs(z^(-(b - a + 2)))
            @test abs(CGL.U(b - a + 2, b + 2, z)) >=
                  0.9CU.U_bmap2_bp2 * abs(z^(-(b - a + 2)))

            # U_dz

            @test abs(CGL.U_dz(a, b, z)) <= CU.U_dz_a_b * abs(z^(-a - 1))
            @test abs(CGL.U_dz(a, b, z)) >= 0.9CU.U_dz_a_b * abs(z^(-a - 1))

            @test abs(CGL.U_dz(b - a, b, z)) <= CU.U_dz_bma_b * abs(z^(-(b - a) - 1))
            @test abs(CGL.U_dz(b - a, b, z)) >= 0.9CU.U_dz_bma_b * abs(z^(-(b - a) - 1))

            # U_dz_dz

            @test abs(CGL.U_dz(a, b, z, 2)) <= CU.U_dz_dz_a_b * abs(z^(-a - 2))
            @test abs(CGL.U_dz(a, b, z, 2)) >= 0.9CU.U_dz_dz_a_b * abs(z^(-a - 2))

            @test abs(CGL.U_dz(b - a, b, z, 2)) <= CU.U_dz_dz_bma_b * abs(z^(-(b - a) - 2))
            @test abs(CGL.U_dz(b - a, b, z, 2)) >=
                  0.9CU.U_dz_dz_bma_b * abs(z^(-(b - a) - 2))

            # U_dz_dz_dz

            @test abs(CGL.U_dz(a, b, z, 3)) <= CU.U_dz_dz_dz_a_b * abs(z^(-a - 3))
            @test abs(CGL.U_dz(a, b, z, 3)) >= 0.9CU.U_dz_dz_dz_a_b * abs(z^(-a - 3))

            @test abs(CGL.U_dz(b - a, b, z, 3)) <=
                  CU.U_dz_dz_dz_bma_b * abs(z^(-(b - a) - 3))
            @test abs(CGL.U_dz(b - a, b, z, 3)) >=
                  0.9CU.U_dz_dz_dz_bma_b * abs(z^(-(b - a) - 3))

            # U_da

            @test abs(CGL.U_da(a, b, z)) <= CU.U_da_a_b * abs(log(z) * z^(-a))
            @test abs(CGL.U_da(a, b, z)) >= 0.9CU.U_da_a_b * abs(log(z) * z^(-a))

            @test abs(CGL.U_da(a + 1, b + 1, z)) <=
                  CU.U_da_ap1_bp1 * abs(log(z) * z^(-(a + 1)))
            @test abs(CGL.U_da(a + 1, b + 1, z)) >=
                  0.9CU.U_da_ap1_bp1 * abs(log(z) * z^(-(a + 1)))

            @test abs(CGL.U_da(a + 2, b + 2, z)) <=
                  CU.U_da_ap2_bp2 * abs(log(z) * z^(-(a + 2)))
            @test abs(CGL.U_da(a + 2, b + 2, z)) >=
                  0.9CU.U_da_ap2_bp2 * abs(log(z) * z^(-(a + 2)))

            @test abs(CGL.U_da(b - a, b, z)) <= CU.U_da_bma_b * abs(log(z) * z^(-(b - a)))
            @test abs(CGL.U_da(b - a, b, z)) >=
                  0.9CU.U_da_bma_b * abs(log(z) * z^(-(b - a)))

            @test abs(CGL.U_da(b - a + 1, b + 1, z)) <=
                  CU.U_da_bmap1_bp1 * abs(log(z) * z^(-(b - a + 1)))
            @test abs(CGL.U_da(b - a + 1, b + 1, z)) >=
                  0.9CU.U_da_bmap1_bp1 * abs(log(z) * z^(-(b - a + 1)))
        end
    end

    @testset "Lemma U_da asymptotic expansion" begin
        # Test some of the formulas used in the lemma. We only make
        # approximate tests. This doesn't prove anything, it is only
        # to reduce the risk of typos.

        κ, ϵ, λ = (Arb(0.917383), Arb(0.01), CGLParams{Arb}(3, 1.0, 1.0, 0.02))
        a, b, c = CGL._abc(κ, ϵ, λ)
        z = c * Arb(10)^2

        ### Check the integral representation of U(a, b, z) with γ
        integrand_U =
            (t; analytic) ->
                exp(-z * t) *
                Arblib.pow_analytic!(zero(t), t, a - 1, analytic) *
                Arblib.pow_analytic!(zero(t), 1 + t, -(a - b + 1), analytic)

        γ = Acb(Arblib.arg!(Arb(), z))

        @test Arblib.integrate(
            integrand_U,
            1e-15,
            100 * exp(-im * γ),
            check_analytic = true,
        ) / gamma(a) ≈ CGL.U(a, b, z) rtol = 1e-6

        ### Check the integral representation for χ

        # Definition of χ
        χ_1(t, n) =
            (1 + t)^(-(a - b + 1)) -
            sum(k -> (-1)^k * rising(a - b + 1, k) / factorial(k) * t^k, 0:n-1)

        # Integral representation of χ
        integrand_χ_2 =
            (t, n) ->
                (w; analytic) ->
                    Arblib.pow_analytic!(zero(w), w - 1, Acb(-(a - b + 1)), analytic) /
                    (w^n * (w + t))
        χ_2(t, n) =
            (-1)^n * sinpi(a - b + 1) / π *
            t^n *
            Arblib.integrate(integrand_χ_2(t, n), 1 + 1e-15, 100, check_analytic = true)

        @test χ_1(10 * exp(-im * γ), 5) ≈ χ_2(10 * exp(-im * γ), 5) rtol = 1e-6

        ### Uniform bound of inv(abs(1 + t / w))
        ρ_γ = γ -> ifelse(abs(γ) < Arb(π) / 2, one(γ), inv(sin(π - abs(γ))))

        for γ in range(Arb(-3), 3, 10)
            for t_div_w in range(0, 10, 1000) * exp(-Acb(0, γ))
                @test inv(abs(1 + t_div_w)) <= ρ_γ(γ)
            end
        end
    end
end
