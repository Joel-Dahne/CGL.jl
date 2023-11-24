@testset "abspow" begin
    x = ArbSeries((1, 2, 3))
    y = Arb(2.3)

    # Check that the rules for the first three terms used in the
    # implementation of CGL.abspow! are correct. To have something to
    # compare it to we use a non-zero x[0].
    res = CGL.abspow(x, y)
    @test Arblib.overlaps(res[0], CGL.abspow(x[0], y))
    @test Arblib.overlaps(res[1], y * CGL.abspow(x[0], y - 1) * x[1])
    @test Arblib.overlaps(
        res[2],
        y * ((y - 1) * CGL.abspow(x[0], y - 2) * x[1]^2 + CGL.abspow(x[0], y - 1) * 2x[2]) /
        2,
    )

    @test Arblib.overlaps(CGL.abspow(x, y), CGL.abspow(x, add_error(y, Mag(1e-15))))
end
