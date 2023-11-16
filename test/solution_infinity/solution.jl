@testset "solution_infinity" begin
    params = [CGL._params.(Arb, 1, d) for d in [1, 3]]

    @testset "solution_infinity" begin
        @testset "Parameters $i" for (i, (μ, γ, κ, ξ₁, λ)) in enumerate(params)
            res_F64_1 = CGL.solution_infinity(
                Complex{Float64}(γ),
                Float64(κ),
                Float64(ξ₁),
                gl_params(Float64, λ),
                order = 1,
            )
            res_F64_2 = CGL.solution_infinity(
                Complex{Float64}(γ),
                Float64(κ),
                Float64(ξ₁),
                gl_params(Float64, λ),
                order = 2,
            )
            res_Arb_1 = CGL.solution_infinity(γ, κ, ξ₁, λ, order = 1)
            res_Arb_2 = CGL.solution_infinity(γ, κ, ξ₁, λ, order = 2)

            @test res_F64_1 ≈ ComplexF64.(res_Arb_1)
            # IMPROVE: Once the second order bounds are improved the
            # rtol should tightened
            @test res_F64_2 ≈ ComplexF64.(res_Arb_2) rtol = 1e-4

            @test all(Arblib.overlaps.(res_Arb_1, res_Arb_2))
        end
    end

    @testset "solution_infinity_jacobian" begin
        @testset "Parameters $i" for (i, (μ, γ, κ, ξ₁, λ)) in enumerate(params)
            # F64 currently only implements first order bounds
            res_jacobian_F64_1 = CGL.solution_infinity_jacobian(
                Complex{Float64}(γ),
                Float64(κ),
                Float64(ξ₁),
                gl_params(Float64, λ),
            )

            res_jacobian_Arb_1 = CGL.solution_infinity_jacobian(γ, κ, ξ₁, λ, order = 1)
            res_jacobian_Arb_2 = CGL.solution_infinity_jacobian(γ, κ, ξ₁, λ, order = 2)

            @test res_jacobian_F64_1 ≈ ComplexF64.(res_jacobian_Arb_1)

            @test all(Arblib.overlaps.(res_jacobian_Arb_1, res_jacobian_Arb_2))
        end
    end
end
