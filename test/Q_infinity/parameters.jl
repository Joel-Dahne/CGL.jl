@testset "parameter" begin
    ξ = Arb(30)
    κ = Arb(0.493223)
    ϵ = Arb(0.1)
    λ = CGLParams{Arb}(1, 1.0, 2.3, 0.2)

    ξF64 = Float64(ξ)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    λF64 = CGLParams{Float64}(λ)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)

    @testset "a" begin
        @test CGL._a(κ, ϵ, λ) ≈ CGL._a(κF64, ϵF64, λF64)

        @test CGL._a_dκ(κ, ϵ, λ) ≈ CGL._a_dκ(κF64, ϵF64, λF64)
        @test Arblib.overlaps(CGL._a_dκ(κ, ϵ, λ), CGL._a(ArbSeries((κ, 1)), ϵ, λ)[1])
        @test CGL._a_dκ(κF64, ϵF64, λF64) ≈ fdm(κ -> CGL._a(κ, ϵF64, λF64), κF64)
    end

    @testset "b" begin
        @test CGL._b(κ, ϵ, λ) ≈ CGL._b(κF64, ϵF64, λF64)
    end

    @testset "c" begin
        @test CGL._c(κ, ϵ, λ) ≈ CGL._c(κF64, ϵF64, λF64)

        @test CGL._c_dκ(κ, ϵ, λ) ≈ CGL._c_dκ(κF64, ϵF64, λF64)
        @test Arblib.overlaps(CGL._c_dκ(κ, ϵ, λ), CGL._c(ArbSeries((κ, 1)), ϵ, λ)[1])
        @test CGL._c_dκ(κF64, ϵF64, λF64) ≈ fdm(κ -> CGL._c(κ, ϵF64, λF64), κF64)

        @test CGL._c_dϵ(κ, ϵ, λ) ≈ CGL._c_dϵ(κF64, ϵF64, λF64)
        @test Arblib.overlaps(CGL._c_dϵ(κ, ϵ, λ), CGL._c(κ, ArbSeries((ϵ, 1)), λ)[1])
        @test CGL._c_dϵ(κF64, ϵF64, λF64) ≈ fdm(ϵ -> CGL._c(κF64, ϵ, λF64), ϵF64)
    end

    @testset "B_W" begin
        @test CGL.B_W(κ, ϵ, λ) ≈ CGL.B_W(κF64, ϵF64, λF64)

        @test CGL.B_W_dκ(κ, ϵ, λ) ≈ CGL.B_W_dκ(κF64, ϵF64, λF64)
        @test Arblib.overlaps(CGL.B_W_dκ(κ, ϵ, λ), CGL.B_W(ArbSeries((κ, 1)), ϵ, λ)[1])
        @test CGL.B_W_dκ(κF64, ϵF64, λF64) ≈ fdm(κ -> CGL.B_W(κ, ϵF64, λF64), κF64)

        @test CGL.B_W_dϵ(κ, ϵ, λ) ≈ CGL.B_W_dϵ(κF64, ϵF64, λF64)
        @test Arblib.overlaps(CGL.B_W_dϵ(κ, ϵ, λ), CGL.B_W(κ, ArbSeries((ϵ, 1)), λ)[1])
        @test CGL.B_W_dϵ(κF64, ϵF64, λF64) ≈ fdm(ϵ -> CGL.B_W(κF64, ϵ, λF64), ϵF64)
    end
end
