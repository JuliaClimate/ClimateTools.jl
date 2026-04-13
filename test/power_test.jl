using ClimateTools
using Test
using YAXArrays
using DimensionalData
using Dates

@testset "Lomb-Scargle" begin
    times = collect(Date(2020, 1, 1):Day(1):Date(2020, 1, 20))
    signal = sin.((collect(1:20) .* 2 .* pi) ./ 5)
    cube = YAXArray((Dim{:time}(times),), signal)

    result = ClimateTools.LombScargle.lombscargle(cube)
    result_array = Array(result)

    date_offsets = getproperty.(Date.(times) .- Date(first(times)), :value)
    plan = ClimateTools.LombScargle.plan(date_offsets, signal)
    periodogram = ClimateTools.LombScargle.lombscargle(plan)

    @test size(result) == (3,)
    @test Tuple(Symbol.(name.(result.axes))) == (:LombScargle,)
    @test result_array[1] ≈ Float64(ClimateTools.LombScargle.M(periodogram))
    @test result_array[2] ≈ Float64(ClimateTools.LombScargle.findmaxperiod(periodogram)[1])
    @test result_array[3] ≈ Float64(ClimateTools.LombScargle.findmaxpower(periodogram))

    sparse_values = Union{Missing, Float64}[1, missing, missing, 2, missing, 3, missing, missing, 4, missing, missing, 5]
    sparse_times = collect(Date(2020, 2, 1):Day(1):Date(2020, 2, 12))
    sparse_cube = YAXArray((Dim{:time}(sparse_times),), sparse_values)
    sparse_result = ClimateTools.LombScargle.lombscargle(sparse_cube)

    @test all(ismissing, Array(sparse_result))
end