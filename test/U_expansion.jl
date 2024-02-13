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

        CU = CGL.C_U(a, b, z₁)
        CU_dz = CGL.C_U_dz(a, b, z₁)
        CU_da = CGL.C_U_da(a, b, z₁)

        for z in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64] .* z₁
            @test abs(CGL.U(a, b, z)) <= CU * abs(z^(-a))
            @test abs(CGL.U(a, b, z)) >= 0.9CU * abs(z^(-a))

            @test abs(CGL.U_dz(a, b, z)) <= CU_dz * abs(z^(-a - 1))
            @test abs(CGL.U_dz(a, b, z)) >= 0.9CU_dz * abs(z^(-a - 1))

            @test abs(CGL.U_da(a, b, z)) <= CU_da * abs(log(z) * z^(-a))
            @test abs(CGL.U_da(a, b, z)) >= 0.9CU_da * abs(log(z) * z^(-a))
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
