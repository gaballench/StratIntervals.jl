using StratIntervals
using Test
using Distributions
using QuadGK
using Random
using Optim

Random.seed!(1993)

@testset "StratIntervals.jl" begin
    @testset "Conflation" begin
        # test conflation code against the theoretical expectation of &(Normal(0, 1), Normal(0, 1)) = Normal(0, sqrt(1/2))
        xx = -4:0.001:4
        @test all(map(x -> conflate(product_distribution([Normal(), Normal()]), x), xx) .≈ map(x -> pdf(Normal(x, sqrt(1/2)), 0.0), xx))
        # test that several arbitrary conflations are PDF because int_{-Inf}^{Inf} &(f(x)_1^N) ≈ 1.0
        # test conflation of Vector{ContinuousUnivariateDistribution}
        integral_normals, err1 = quadgk(x -> conflate([Normal(0,1), Normal(3,2), Normal(5,4)], x), -Inf, Inf, rtol=1e-8)
        @test integral_normals ≈ 1.0
        # test conflation of Array{Float64, 2}
        distribs_matrix = [rand(Exponential(5), 50) rand(Uniform(2, 10), 50)]
        integral_empiricals, err2 = quadgk(x -> conflate(distribs_matrix, x), -Inf, Inf, rtol=1e-8)
        @test integral_empiricals ≈ 1.0
        # test that quant(taus) works
        normals_vec1 = rand(Normal(0.0, 1.0), 500000)
        normals_vec2 = rand(Normal(0.0, 1.0), 500000)
        normals_matrix = [normals_vec1 normals_vec1]
        @test isapprox(quant(normals_matrix, 0.5, true, -10.0, 10.0)[1], 0.0, atol=1e-2)
        # testing quant(taus) with a single distribution in a 1-column matrix with different theoretical quantiles
        @test isapprox(quant(normals_matrix[:, 1:1], 0.5, true, -10.0, 10.0)[1], median(Normal(0, 1)), atol=1e-2)
        @test isapprox(quant(normals_matrix[:, 1:1], 0.25, true, -10.0, 10.0)[1], quantile(Normal(0, 1), 0.25), atol=1e-2)
        @test isapprox(quant(normals_matrix[:, 1:1], 0.75, true, -10.0, 10.0)[1], quantile(Normal(0, 1), 0.75), atol=1e-2)
    end
    @testset "Distributions" begin
        # tests here
    end
    @testset "Stratigraphic intervals" begin
        # tests here
    end
end
