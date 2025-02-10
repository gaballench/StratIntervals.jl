using StratIntervals
using Test
using Distributions
using QuadGK
using Random

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
    end
    @testset "Distributions" begin
        # tests here
    end
    @testset "Stratigraphic intervals" begin
        # tests here
    end
end
