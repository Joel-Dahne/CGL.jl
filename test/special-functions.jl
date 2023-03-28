@testset "special-functions" begin
    @testset "rising" begin
        @test rising(Arb(5), Arb(3)) == 210
        @test rising(Arb(5), 3) == 210
        @test rising(ArbSeries((5, 1)), 3) == ArbSeries((210, 107))
        @test rising(5.0, 3.0) == 210
        @test rising(5.0, 3) == 210
        @test rising(5, 3) == 210
    end

    @testset "hypgeom_u" begin
        @test Float64(hypgeom_u(Arb(0.1), Arb(0.2), Arb(0.3))) == hypgeom_u(0.1, 0.2, 0.3)
        @test ComplexF64(hypgeom_u(Acb(0.1, 1.1), Acb(0.2, 1.2), Acb(0.3, 1.3))) ==
              hypgeom_u(0.1 + im * 1.1, 0.2 + im * 1.2, 0.3 + im * 1.3)

        # For ArbSeries it is convenient to test it by checking that
        # it gives correct results by expanding in a Taylor series and
        # evaluating
        let a = Arb(0.1), b = Arb(0.2)
            for z0 in Arb[0.2, 0.5, 1.5]
                for r in Mag[1e-10, 1e-5, 1e-1]
                    z = add_error(z0, r)

                    # Test hypgeom_u(a, b, z)
                    res1 = hypgeom_u(a, b, z0 - r)
                    res2 = hypgeom_u(a, b, z0 + r)
                    for degree = 0:5
                        # Compute Taylor expansion and remainder term
                        p = hypgeom_u(a, b, ArbSeries((z0, 1); degree))
                        remainder = ArbExtras.taylor_remainder(
                            hypgeom_u(a, b, ArbSeries((z, 1), degree = degree + 1)),
                            z,
                        )

                        @test Arblib.overlaps(p(-Arb(r)) + remainder, res1)
                        @test Arblib.overlaps(p(Arb(r)) + remainder, res2)
                    end

                    # Test hypgeom_u(a, b, atan(z))
                    res1 = hypgeom_u(a, b, atan(z0 - r))
                    res2 = hypgeom_u(a, b, atan(z0 + r))
                    for degree = 0:5
                        p = hypgeom_u(a, b, atan(ArbSeries((z0, 1); degree)))
                        remainder = ArbExtras.taylor_remainder(
                            hypgeom_u(a, b, atan(ArbSeries((z, 1), degree = degree + 1))),
                            z,
                        )
                        @test Arblib.overlaps(p(-Arb(r)) + remainder, res1)
                        @test Arblib.overlaps(p(Arb(r)) + remainder, res2)
                    end
                end
            end
        end
    end
end
