using ClimateTools
using Test
using NetCDF
using YAXArrays
using DimensionalData
using Statistics
using Dates

@testset "ClimateTools.jl" begin    
    # Write your tests here.    

    obsfile = joinpath(dirname(@__FILE__), "data", "obs.nc")
    reffile = joinpath(dirname(@__FILE__), "data", "ref.nc")
    futfile = joinpath(dirname(@__FILE__), "data", "fut.nc")

    obs = Cube(YAXArrays.open_dataset(obsfile))
    ref = Cube(YAXArrays.open_dataset(reffile))
    fut = Cube(YAXArrays.open_dataset(futfile))

   

    qq = qqmap(obs, ref, fut, detrend=true, method="multiplicative")
    @test round(qq[45,3,3], digits=5) == -7.85231
    @test round(mean(qq), digits=5) == 0.39328
    @test round(std(qq), digits=5) == 6.3777
    
    

    qq = qqmap(obs, ref, fut, detrend=true, method="additive")
    @test round(qq[45,3,3], digits=5) == -7.83775
    @test round(mean(qq), digits=5) == 0.25184
    @test round(std(qq), digits=5) == 5.54588
    
    # detrend false, multiplicative
    qq = qqmap(obs, ref, fut, detrend=false, method="multiplicative")      
    @test round(qq[45,3,3], digits=5) == -7.82242
    @test round(mean(qq), digits=5) == 0.38062
    @test round(std(qq), digits=5) == 6.36466
    
    # detrend false, additive
    qq = qqmap(obs, ref, fut, detrend=false, method="additive")
    @test round(qq[45,3,3], digits=5) == -7.80781
    @test round(mean(qq), digits=5) == 0.23945
    @test round(std(qq), digits=5) == 5.51853

    # Plots.density(obs[:,3,3].data, linewidth=3,label="OBS")
    # Plots.density!(ref[:,3,3].data, linewidth=2.5,label="REF")
    # Plots.density!(fut[:,3,3].data, linewidth=2.5,label="FUT")
    # Plots.density!(qq[:,3,3].data, linewidth=2.5,label="QQ")

end
