@testset "verify_and_refine_root" begin
    f = x -> [sin(x[1]), sin(x[2]), sin(x[3]), cos(x[4])]
    df = x -> [cos(x[1]) 0 0 0; 0 cos(x[2]) 0 0; 0 0 cos(x[3]) 0; 0 0 0 -sin(x[4])]

    root_true = SVector(Arb(0), Arb(π), 2Arb(π), Arb(π) / 2)
    root_enclosure = add_error.(root_true, Mag(1))

    root = CGL.verify_and_refine_root(f, df, root_enclosure)

    atol = 0
    rtol = 4eps(one(first(root)))

    @test all(Arblib.overlaps.(root, root_true))
    @test all(ArbExtras.check_tolerance.(root; atol, rtol))
end
