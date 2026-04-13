using ClimateTools
using Test
using YAXArrays
using DimensionalData
using Dates

@testset "Climatology" begin
    times = collect(Date(2000, 1, 1):Month(1):Date(2001, 12, 1))
    data = reshape(Float64.(repeat(1:12, 2)), 1, 1, :)
    cube = YAXArray((Dim{:longitude}([0.0]), Dim{:latitude}([45.0]), Dim{:time}(times)), data)

    result = climato_tp(cube; fct=sum, iduree=3, lead=0)
    result_array = Array(result)

    @test size(result) == (2, 1, 1, 12)
    @test Tuple(Symbol.(name.(result.axes))) == (:time, :longitude, :latitude, :Mois)
    @test collect(lookup(result, :Mois)) == 1:12
    @test year.(collect(lookup(result, :time))) == [2000, 2001]

    @test result_array[:, 1, 1, 1] == [6.0, 6.0]
    @test result_array[:, 1, 1, 11] == [24.0, 24.0]
    @test result_array[:, 1, 1, 12] == [15.0, 15.0]
end