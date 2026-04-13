using ClimateTools
using Test
using YAXArrays
using DimensionalData
using Dates

@testset "Markov Switching" begin
    times = collect(Date(2020, 1, 1):Day(1):Date(2020, 2, 29))
    values = vcat(fill(-2.0, 30), fill(2.0, 30))
    cube = YAXArray((Dim{:time}(times),), values)

    result = ClimateTools.MSModel(cube)
    @test size(result) == (1,)
    @test Tuple(Symbol.(name.(result.axes))) == (:MSM,)
    @test Array(result)[1] isa ClimateTools.MarSwitching.MSM

    result_three = ClimateTools.MSModel(cube; k_regime=3)
    @test size(result_three) == (1,)
    @test Array(result_three)[1] isa ClimateTools.MarSwitching.MSM
end