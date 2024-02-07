@testset "special-functions" begin
    @testset "rising" begin
        @test rising(Arb(5), Arb(3)) == 210
        @test rising(Arb(5), 3) == 210
        @test rising(ArbSeries((5, 1)), 3) == ArbSeries((210, 107))
        @test rising(5.0, 3.0) == 210
        @test rising(5.0, 3) == 210
        @test rising(5, 3) == 210
    end

    @testset "U" begin
        @test Float64(U(Arb(0.1), Arb(0.2), Arb(0.3))) == U(0.1, 0.2, 0.3)
        @test ComplexF64(U(Acb(0.1, 1.1), Acb(0.2, 1.2), Acb(0.3, 1.3))) ==
              U(0.1 + im * 1.1, 0.2 + im * 1.2, 0.3 + im * 1.3)

        # For ArbSeries it is convenient to test it by checking that
        # it gives correct results by expanding in a Taylor series and
        # evaluating
        let a = Arb(0.1), b = Arb(0.2)
            for z0 in Arb[0.2, 0.5, 1.5]
                for r in Mag[1e-10, 1e-5, 1e-1]
                    z = add_error(z0, r)

                    # Test U(a, b, z)
                    res1 = U(a, b, z0 - r)
                    res2 = U(a, b, z0 + r)
                    for degree = 0:5
                        # Compute Taylor expansion and remainder term
                        p = U(a, b, ArbSeries((z0, 1); degree))
                        remainder = ArbExtras.taylor_remainder(
                            U(a, b, ArbSeries((z, 1), degree = degree + 1)),
                            z,
                        )

                        @test Arblib.overlaps(p(-Arb(r)) + remainder, res1)
                        @test Arblib.overlaps(p(Arb(r)) + remainder, res2)
                    end

                    # Test U(a, b, atan(z))
                    res1 = U(a, b, atan(z0 - r))
                    res2 = U(a, b, atan(z0 + r))
                    for degree = 0:5
                        p = U(a, b, atan(ArbSeries((z0, 1); degree)))
                        remainder = ArbExtras.taylor_remainder(
                            U(a, b, atan(ArbSeries((z, 1), degree = degree + 1))),
                            z,
                        )
                        @test Arblib.overlaps(p(-Arb(r)) + remainder, res1)
                        @test Arblib.overlaps(p(Arb(r)) + remainder, res2)
                    end
                end
            end
        end
    end

    ξ = Arb(30)
    κ = Arb(0.493223)
    ϵ = Arb(0.1)
    λ = CGLParams{Arb}(1, 1.0, 2.3, 0.2)
    a, b, c = CGL._abc(κ, ϵ, λ)
    z = c * ξ^2

    aF64 = ComplexF64(a)
    bF64 = ComplexF64(b)
    zF64 = ComplexF64(z)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)
    fdm2 = central_fdm(5, 2)

    @testset "U_dz" begin
        @test U_dz(a, b, z) ≈ fdm(z_r -> U(aF64, bF64, z_r + im * imag(zF64)), real(zF64)) rtol =
            1e-12

        @test U_dz(a, b, z, 2) ≈
              fdm2(z_r -> U(aF64, bF64, z_r + im * imag(zF64)), real(zF64)) rtol = 1e-7
    end

    @testset "U_da" begin
        @test U_da(a, b, z) ≈ fdm(a_r -> U(a_r + im * imag(aF64), bF64, zF64), real(aF64)) rtol =
            1e-12


        # Compare with U_1f1. This only works well for exact
        # input at very high precision though.
        a_s = AcbSeries((midpoint(a), 1))
        b_s = AcbSeries(midpoint(b), degree = 1)
        z_s = AcbSeries(midpoint(z), degree = 1)
        res = zero(a_s)

        # This version only works well at very high precision and
        # not for wide input.
        Arblib.hypgeom_u_1f1_series!(res, a_s, b_s, z_s, 2, 10precision(a_s))
        @test Arblib.overlaps(
            res[1],
            U_da(midpoint(Acb, a), midpoint(Acb, b), midpoint(Acb, z)),
        )

        Arblib.hypgeom_u_1f1_series!(res, a_s, b_s, 2z_s, 2, 10precision(a_s))
        @test Arblib.overlaps(
            res[1],
            U_da(midpoint(Acb, a), midpoint(Acb, b), midpoint(Acb, 2z)),
        )
    end

    @testset "U_dzda" begin
        @test U_dzda(a, b, z) ≈ fdm(
            a_r -> fdm(
                z_r -> U(a_r + im * imag(aF64), bF64, z_r + im * imag(zF64)),
                real(zF64),
            ),
            real(aF64),
        ) rtol = 1e-9
        @test U_dzda(a, b, z) ≈ fdm(
            z_r -> fdm(
                a_r -> U(a_r + im * imag(aF64), bF64, z_r + im * imag(zF64)),
                real(aF64),
            ),
            real(zF64),
        ) rtol = 1e-8

        # Here we get very large errors for the fdm. But it at least
        # checks the first few digits.
        @test U_dzda(a, b, z, 2) ≈ fdm(
            a_r -> fdm2(
                z_r -> U(a_r + im * imag(aF64), bF64, z_r + im * imag(zF64)),
                real(zF64),
            ),
            real(aF64),
        ) rtol = 1e-3
        @test U_dzda(a, b, z, 2) ≈ fdm2(
            z_r -> fdm(
                a_r -> U(a_r + im * imag(aF64), bF64, z_r + im * imag(zF64)),
                real(aF64),
            ),
            real(zF64),
        ) rtol = 1e-2
    end
end
