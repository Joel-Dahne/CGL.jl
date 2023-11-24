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
            res_Arb = CGL.solution_infinity(γ, κ, ξ₁, λ)

            @test res_F64_1 ≈ ComplexF64.(res_Arb) rtol = 1e-4
            @test res_F64_2 ≈ ComplexF64.(res_Arb) rtol = 1e-6
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

            res_jacobian_Arb = CGL.solution_infinity_jacobian(γ, κ, ξ₁, λ)

            # For the first column we get good enclosures
            @test res_jacobian_F64_1[:, 1] ≈ ComplexF64.(res_jacobian_Arb[:, 1])

            # For the second column we get very bad enclosures and
            # hence not very good agreement. We only check that the
            # F64 approximation is inside the compute enclosure.
            @test all(
                Arblib.contains.(res_jacobian_Arb[:, 2], Acb.(res_jacobian_F64_1[:, 2])),
            )
        end
    end
end
