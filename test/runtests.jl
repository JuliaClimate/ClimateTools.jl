using Climat
using Test
using NetCDF
using YAXArrays
using DimensionalData
# using StableRNGs
using Statistics
using Dates

@testset "Climat.jl" begin    
    # Write your tests here.

    # obsfile = joinpath("/gpfs/home/dl2594/.julia/dev/Climat/test/", "data", "obs.nc")
    # reffile = joinpath("/gpfs/home/dl2594/.julia/dev/Climat/test/", "data", "ref.nc")
    # futfile = joinpath("/gpfs/home/dl2594/.julia/dev/Climat/test/", "data", "fut.nc")

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

# # Generation of random data
# datevec = DateTime(1981,1,1):Dates.Day(1):DateTime(2010,12,31)
# datevec2 = DateTime(2011,1,1):Dates.Day(1):DateTime(2040,12,31)

# lon = [-75.0, -72.5, -70.0, -67.5]
# lat = [35.0, 37.5, 40.0, 42.5]

# mu_obs = 0.0;sigma_obs = 5.0
# mu_ref = 3.15;sigma_ref = 7.4
# mu_fut = 3.5;sigma_fut = 8.0

# obsdata = rand(Normal(mu_obs, sigma_obs), length(datevec), 4, 4)
# refdata = rand(Normal(mu_ref, sigma_ref), length(datevec), 4, 4)
# futdata = rand(Normal(mu_fut, sigma_fut), length(datevec2), 4, 4)


# axlist = (
#     Dim{:Ti}(datevec),
#     Dim{:lon}(lon),
#     Dim{:lat}(lat),
# )

# axlist2 = (
#     Dim{:Ti}(datevec2),
#     Dim{:lon}(lon),
#     Dim{:lat}(lat),
# )

# properties_obs = Dict("mu" => mu_obs, "sigma" => sigma_obs)
# obs = YAXArray(axlist, obsdata, properties_obs)

# properties_ref = Dict("mu" => mu_ref, "sigma" => sigma_ref)
# ref = YAXArray(axlist, refdata, properties_ref)

# properties_fut = Dict("mu" => mu_fut, "sigma" => sigma_fut)
# fut = YAXArray(axlist2, futdata, properties_fut)
