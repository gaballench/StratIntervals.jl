using StratIntervals
using Test
using Distributions

@testset "StratIntervals.jl" begin
    @testset "Conflation" begin
        # test conflation code against the theoretical expectation of &(Normal(0, 1), Normal(0, 1)) = Normal(0, sqrt(1/2))
        xx = -4:0.001:4
        @test all(map(x -> conflate(product_distribution([Normal(), Normal()]), x), xx) .≈ map(x -> pdf(Normal(x, sqrt(1/2)), 0.0), xx))
    end
    @testset "Distributions" begin
        # tests here
    end
    @testset "Stratigraphic intervals" begin
        # tests here
    end
end
