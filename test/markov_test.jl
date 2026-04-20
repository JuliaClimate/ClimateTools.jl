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

@testset "Markov Switching multidimensional cube" begin
    times = collect(Date(2021, 1, 1):Day(1):Date(2021, 2, 19))
    lon = [10.0, 20.0]
    lat = [45.0, 50.0]

    base_signal = vcat(fill(-3.0, 25), fill(3.0, 25))
    data = Array{Float64}(undef, length(lon), length(lat), length(times))
    for ilon in eachindex(lon), ilat in eachindex(lat)
        data[ilon, ilat, :] .= base_signal .+ 0.1 * ilon .- 0.05 * ilat
    end

    cube = YAXArray((Dim{:longitude}(lon), Dim{:latitude}(lat), Dim{:time}(times)), data)
    result = ClimateTools.MSModel(cube; k_regime=2, intercept="switching")

    @test size(result) == (length(lon), length(lat), 1)
    @test Tuple(Symbol.(name.(result.axes))) == (:longitude, :latitude, :MSM)
    @test all(model -> model isa ClimateTools.MarSwitching.MSM, Array(result))
end

@testset "Markov Switching validation" begin
    cube = YAXArray((Dim{:longitude}([1.0, 2.0]), Dim{:latitude}([3.0, 4.0])), [1.0 2.0; 3.0 4.0])
    @test_throws ErrorException ClimateTools.MSModel(cube)
end