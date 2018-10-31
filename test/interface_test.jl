using Unitful: K, °C, m, mm, s, kg
using Unitful: @u_str, ustrip, uconvert

@testset "Interface" begin

# test that load return a ClimGrid type
file1 = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
file2 = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
files = [file1, file2]
C = load(files, "tas")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
C = load(filenc, "tas")

@test load(filenc, "tas", data_units = "Celsius")[2] == "Celsius"
@test load(filenc, "pr", data_units = "mm")[2] == "mm"
@test typeof(C) == ClimGrid{AxisArrays.AxisArray{Unitful.Quantity{Float32,Unitful.Dimensions{(Unitful.Dimension{:Temperature}(1//1),)},Unitful.FreeUnits{(Unitful.Unit{:Kelvin,Unitful.Dimensions{(Unitful.Dimension{:Temperature}(1//1),)}}(0, 1//1),),Unitful.Dimensions{(Unitful.Dimension{:Temperature}(1//1),)}}},3,Array{Unitful.Quantity{Float32,Unitful.Dimensions{(Unitful.Dimension{:Temperature}(1//1),)},Unitful.FreeUnits{(Unitful.Unit{:Kelvin,Unitful.Dimensions{(Unitful.Dimension{:Temperature}(1//1),)}}(0, 1//1),),Unitful.Dimensions{(Unitful.Dimension{:Temperature}(1//1),)}}},3},Tuple{AxisArrays.Axis{:lon,Array{Float32,1}},AxisArrays.Axis{:lat,Array{Float32,1}},AxisArrays.Axis{:time,Array{Dates.DateTime,1}}}}}


fileorog = joinpath(dirname(@__FILE__), "data", "orog_fx_GFDL-ESM2G_historicalMisc_r0i0p0.nc")
orog = load2D(fileorog, "orog")
@test size(orog[1]) == (144, 90)
@test orog.typeofvar == "orog"
status, figh = mapclimgrid(orog); @test status == true; PyPlot.close();

filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
polyshp = read(filename,Shapefile.Handle)
P = shapefile_coords_poly(polyshp.shapes[1])
orog = load2D(fileorog, "orog", poly = P)
@test size(orog[1]) == (13, 8)
@test orog.typeofvar == "orog"

# INTERFACE
# B = vcat(C, C)
# @test size(B.data) == (256, 128, 2) # vcat does not look at dimensions
B = merge(C, C)
@test size(B.data) == (256, 128, 1) # C being similar, they should not add up, as opposed to vcat
# Operators +, -, *, /
B = C + C; @test B[1].data[1, 1, 1] == 438.4457f0u"K"
B = C * C; @test B[1].data[1, 1, 1] == 48058.66f0u"K^2"
B = C / C; @test B[1].data[1, 1, 1] == 1.0f0
B = C - C; @test B[1].data[1, 1, 1] == 0.0f0u"K"
B = C - 1.0u"K"; @test B[1].data[1, 1, 1] == 218.2228546142578u"K"
B = C - 1u"K"; @test B[1].data[1, 1, 1] == 218.22285f0u"K"
B = C / 2; @test B[1].data[1, 1, 1] == 109.61143f0u"K"
B = C / 2u"K"; @test B[1].data[1, 1, 1] == 109.61143f0
B = C / 2.2; @test B[1].data[1, 1, 1] == 99.6467520973899u"K"
B = C * 2; @test B[1].data[1, 1, 1] == 438.4457f0u"K"
B = C * 2u"K"; @test B[1].data[1, 1, 1] == 438.4457f0u"K^2"
B = C * 2.2; @test B[1].data[1, 1, 1] == 482.2902801513672u"K"

@test ClimateTools.mean(C) == 278.6421f0u"K"
@test ClimateTools.maximum(C) == 309.09613f0u"K"
@test ClimateTools.minimum(C) == 205.24321f0u"K"
@test ClimateTools.std(C) == 21.92836f0u"K"
@test round(ustrip(ClimateTools.var(C)), digits=3) == 480.853f0

# @test typeof(show(C)) == Dict{Any, Any}
# @test typeof(C[1].data) == Array{Unitful.Quantity{Float32,Unitful.Dimensions{(Unitful.Dimension{:Temperature}(1//1),)},Unitful.FreeUnits{(Unitful.Unit{:Kelvin,Unitful.Dimensions{(Unitful.Dimension{:Temperature}(1//1),)}}(0, 1//1),),Unitful.Dimensions{(Unitful.Dimension{:Temperature}(1//1),)}}},3}

@test C[2] == "K"
@test C[3] == "N/A"
@test C[4] == "720 ppm stabilization experiment (SRESA1B)"
@test C[5] == "N/A"
@test C[6] == "degrees_east"
@test C[7] == "degrees_north"
@test C[8] == joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
@test C[9] == "tas"
@test annualmax(C)[9] == "annualmax"
@test C[10] == "tas"
@test C[11] == "noleap"
@test typeof(C[12]) == Dict{Any, Any}
@test C[12]["project_id"] == "IPCC Fourth Assessment"
@test_throws ErrorException C[13]
@test annualmax(C)[10] == "tas"
@test size(C) == (22, )
@test size(C, 1) == 22
@test length(C) == 22
# @test endof(C) == 21
@test_throws MethodError C[end]
@test ndims(C) == 1


# Spatial subset
filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
polyshp = read(filename,Shapefile.Handle)
x, y = shapefile_coords(polyshp.shapes[1])
P = [x y]
P = P'
C = load(filenc, "tas")
Csub = spatialsubset(C, P)
@test size(Csub[1]) == (23, 12, 1)
@test Csub[1][1, 1, 1] == 294.6609f0u"K"
Csub = spatialsubset(C, P')
@test size(Csub[1]) == (23, 12, 1)
@test Csub[1][1, 1, 1] == 294.6609f0u"K"
C = load(filenc, "ua")
Csub = spatialsubset(C, P)
@test size(Csub[1]) == (23, 12, 17, 1)
@test Csub[1][12, 1, 1, 1] == 6.658482f0u"m/s"
@test isnan(Csub[1][1, 1, 1, 1])

poly= [[NaN 10 -10 -10 10 10];[NaN -10 -20 10 10 -10]] # meridian test
C = load(filenc, "tas", poly=poly)

# Spatial subset
C = load(filenc, "tas")
Csub = temporalsubset(C, (2000, 05, 16, 12, 0, 0), (2000, 05, 16, 12, 0, 0))
@test Csub[1][1, 1, 1] == 219.22285f0
@test Csub[1][Axis{:time}][1] == DateTimeNoLeap(2000, 05, 16, 12)
B = load(filenc, "tas", start_date=(2000, 05, 16, 12), end_date=(2000, 05, 16, 12))
@test B[1] == C[1]

# Time resolution
timevec = [1, 2, 3]
@test ClimateTools.timeresolution(timevec) == "24h"
@test ClimateTools.daymean_factor(ClimateTools.timeresolution(timevec)) == 1

timevec = [1.0, 1.5, 2.0]
@test ClimateTools.timeresolution(timevec) == "12h"
@test ClimateTools.daymean_factor(ClimateTools.timeresolution(timevec)) == 2

timevec = [1.25, 1.5, 1.75]
@test ClimateTools.timeresolution(timevec) == "6h"
@test ClimateTools.daymean_factor(ClimateTools.timeresolution(timevec)) == 4

timevec = [1.125, 1.25, 1.375]
@test ClimateTools.timeresolution(timevec) == "3h"
@test ClimateTools.daymean_factor(ClimateTools.timeresolution(timevec)) == 8

timevec = [1.0, 1.04167, 1.08334]
@test ClimateTools.timeresolution(timevec) == "1h"
@test ClimateTools.daymean_factor(ClimateTools.timeresolution(timevec)) == 24

timevec = [0.0, 365.0, 365.0*2]
@test ClimateTools.timeresolution(timevec) == "Yearly"

timevec = NetCDF.ncread(filenc, "time")
@test ClimateTools.timeresolution(timevec) == "N/A"
@test ClimateTools.daymean_factor(ClimateTools.timeresolution(timevec)) == 1



# MESHGRID
YV = [1, 2, 3]
XV = [4, 5, 6]
@test meshgrid(XV, YV) == ([4 5 6; 4 5 6; 4 5 6], [1 1 1; 2 2 2; 3 3 3])
@test meshgrid(XV) == ([4 5 6; 4 5 6; 4 5 6], [4 4 4; 5 5 5; 6 6 6])
Q = Array{Int64, 3}(undef, 3, 3, 3)
R = Array{Int64, 3}(undef, 3, 3, 3)
S = Array{Int64, 3}(undef, 3, 3, 3)
Q[:, :, 1] .= [4 5 6; 4 5 6; 4 5 6]
Q[:, :, 2] .= [4 5 6; 4 5 6; 4 5 6]
Q[:, :, 3] .= [4 5 6; 4 5 6; 4 5 6]
R[:, :, 1] .= [1 1 1; 2 2 2; 3 3 3]
R[:, :, 2] .= [1 1 1; 2 2 2; 3 3 3]
R[:, :, 3] .= [1 1 1; 2 2 2; 3 3 3]
S[:, :, 1] .= [4 4 4; 4 4 4; 4 4 4]
S[:, :, 2] .= [5 5 5; 5 5 5; 5 5 5]
S[:, :, 3] .= [6 6 6; 6 6 6; 6 6 6]

@test meshgrid(XV, YV, XV) == (Q, R, S)


@test ndgrid(XV, YV) == ([4 4 4; 5 5 5; 6 6 6], [1 2 3; 1 2 3; 1 2 3])
@test ndgrid([XV, YV, XV]) ==   [[4, 5, 6], [1, 2, 3], [4, 5, 6]]
@test ndgrid(XV) == [4, 5, 6]

# isdefined
C = 1
@test ClimateTools.@isdefined C
@test (ClimateTools.@isdefined T) == false


## INPOLY
@test ClimateTools.leftorright(0.5,0.5, 1,0,1,1) == -1
@test ClimateTools.leftorright(1.5,.5, 1,0,1,1) == 1
@test ClimateTools.leftorright(1,0.5, 1,0,1,1) == 0

poly = Float64[0 0
               0 1
               1 1
               1 0
               0 0]'
p1 = [0.5, 0.5]
p2 = [0.5, 0.99]
p22 = [0.5, 1] # on edge
p23 = [0.5, 0] # on edge
p24 = [0, 0]   # on corner
p25 = [0, .4]   # on edge
p3 = [0.5, 1.1]

@test inpoly(p1, poly)
@test inpoly(p2, poly)
@test inpoly(p22, poly)
@test inpoly(p23, poly)
@test inpoly(p24, poly)
@test inpoly(p25, poly)
@test !inpoly(p3, poly)

# clockwise poly
poly = Float64[0 0
               1 0
               1 1
               0 1
               0 0]'

@test inpoly(p1, poly)
@test inpoly(p2, poly)
@test inpoly(p22, poly)
@test inpoly(p23, poly)
@test inpoly(p24, poly)
@test inpoly(p25, poly)
@test !inpoly(p3, poly)


# cross-over poly
poly = Float64[0 0
               1 0
               0 1
               1 1
               0 0]'
eval(:(@test_broken inpoly(p1, poly) )) # should be true
@test inpoly(p2, poly)
@test inpoly(p22, poly)
@test inpoly(p23, poly)
@test inpoly(p24, poly)
@test !inpoly(p25, poly) # different
@test !inpoly(p3, poly)


# with interior region
poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside interior poly: i.e. labeled as outside
@test !inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               # in-interior
               0.4 0.4
               0.4 0.2
               0.2 0.2
               0.2 0.4
               0.4 0.4
               # interior
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside in-interior poly
@test inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior
               0.1 0.1
               0.1 0.6
               0.6 0.6
               # in-interior
               0.4 0.4
               0.2 0.4
               0.2 0.2
               0.4 0.2
               0.4 0.4
               # interior
               0.6 0.6
               0.6 0.1
               0.1 0.1
               # exterior
               0 0
               0 1
               1 1
               1 0
               0 0]'
# inside in-interior poly
@test inpoly([0.3,0.3], poly)
@test !inpoly([0.3,0.5], poly)

poly = Float64[0 0
               # interior #1
               0.1 0.1
               0.1 0.6
               0.4 0.6
               0.4 0.6
               0.4 0.1
               0.1 0.1
               0 0
               # interior #2
               0.6 0.4
               0.6 0.6
               0.8 0.6
               0.8 0.4
               0.6 0.4
               0 0
               # exterior
               0 1
               1 1
               1 0
               0 0]'
@test !inpoly([0.2,0.4], poly)
@test !inpoly([0.3,0.15], poly)
@test inpoly([0.5,0.4], poly)
@test inpoly([0.5,0.2], poly)
@test !inpoly([0.7,0.5], poly)

# Test Interpolation/regridding
filename = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
C = load(filename, "tas")
# Get lat lon vector
lat = Float32.(C[1][Axis{:lat}][:])
lon = Float32.(C[1][Axis{:lon}][:])
latgrid = Float32.(C.latgrid)
longrid = Float32.(C.longrid)
# Shift longitude by 1
lon .+= Float32(1.0)
longrid .+= Float32(1.0)
axisdata = AxisArray(C[1].data, Axis{:lon}(lon), Axis{:lat}(lat), Axis{:time}(C[1][Axis{:time}][:]))
C2 = ClimGrid(axisdata, variable = "tas", longrid=longrid, latgrid=latgrid, msk=C.msk)
@test regrid(C, C2)[1].data[1, 1, 1] == 219.2400638156467
@test regrid(C, C2, min=0.0, max=0.0)[1].data[1, 1, 1] == 0.0
@test regrid(C, lon, lat)[1].data[1, 1, 1] == 219.2400638156467
@test regrid(C, lon, lat, min=0.0, max=0.0)[1].data[1, 1, 1] == 0.0

# Test applymask
# 1-D data
data = randn(3)
mask = [NaN; 1.;1.]
@test isnan(applymask(data, mask)[1])
@test applymask(data, mask)[2] == data[2]
@test applymask(data, mask)[3] == data[3]

# 2-D data
data = randn(3, 2)
mask = [[NaN; 1; 1] [1.; NaN;1.]]
@test isnan(applymask(data, mask)[1, 1]) && isnan(applymask(data, mask)[2, 2])
@test applymask(data, mask)[2, 1] == data[2, 1]
@test applymask(data, mask)[1, 2] == data[1, 2]
@test applymask(data, mask)[3, 1] == data[3, 1]
@test applymask(data, mask)[3, 2] == data[3, 2]

# 3-D data
data = randn(3, 2, 3)
mask = [[NaN; 1; 1] [1.; NaN;1.]]
@test isnan(applymask(data, mask)[1, 1, 1]) && isnan(applymask(data, mask)[2,2,1]) && isnan(applymask(data, mask)[1, 1, 2]) && isnan(applymask(data, mask)[2, 2, 2]) && isnan(applymask(data, mask)[1,1,3]) && isnan(applymask(data, mask)[2, 2, 3])

for i = 1:size(data, 1)
    @test applymask(data, mask)[2,1,i] == data[2,1,i]
    @test applymask(data, mask)[1,2,i] == data[1,2,i]
    @test applymask(data, mask)[3,1,i] == data[3,1,i]
    @test applymask(data, mask)[3,2,i] == data[3,2,i]
end

# 4-D data
data = randn(3,2,1,3)
mask = [[NaN; 1; 1] [1.; NaN;1.]]
@test isnan(applymask(data, mask)[1, 1, 1, 1]) && isnan(applymask(data, mask)[2,2,1,1]) && isnan(applymask(data, mask)[1, 1, 1,2]) && isnan(applymask(data, mask)[2, 2, 1,2]) && isnan(applymask(data, mask)[1, 1, 1,3]) && isnan(applymask(data, mask)[2, 2, 1,3])

for i = 1:size(data, 1)
    @test applymask(data, mask)[2,1,1,i] == data[2,1,1,i]
    @test applymask(data, mask)[1,2,1,i] == data[1,2,1,i]
    @test applymask(data, mask)[3,1,1,i] == data[3,1,1,i]
    @test applymask(data, mask)[3,2,1,i] == data[3,2,1,i]
end

# Test timeresolution and pr_timefactor
filename = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
timevec = NetCDF.ncread(filename, "time")
@test ClimateTools.pr_timefactor(ClimateTools.timeresolution(timevec)) == 1.0
@test ClimateTools.pr_timefactor("24h") == 86400.0
@test ClimateTools.pr_timefactor("12h") == 43200.0
@test ClimateTools.pr_timefactor("6h") == 21600.0
@test ClimateTools.pr_timefactor("3h") == 10800.0
@test ClimateTools.pr_timefactor("1h") == 3600.0

# timeindex tests
d = collect(DateTime(2000, 01, 01):Day(1):DateTime(2005, 12, 31))
T = typeof(d[1])

start_date = (Inf,)
end_date = (Inf,)
idxtimebeg, idxtimeend = ClimateTools.timeindex(d, start_date, end_date, T)
@test idxtimebeg == 1
@test idxtimeend == length(d)

start_date = (2001, 01, 01)
end_date = (Inf,)
idxtimebeg, idxtimeend = ClimateTools.timeindex(d, start_date, end_date, T)
@test idxtimebeg == 367
@test idxtimeend == length(d)

start_date = (2001, 01, 01)
end_date = (2003, 12, 31)
idxtimebeg, idxtimeend = ClimateTools.timeindex(d, start_date, end_date, T)
@test idxtimebeg == 367
@test idxtimeend == 1461

@test_throws ArgumentError ClimateTools.timeindex(d, end_date, start_date, T)

start_date = (2001,)
end_date = (2004,)
idxtimebeg, idxtimeend = ClimateTools.timeindex(d, start_date, end_date, T)
@test idxtimebeg == 367
@test idxtimeend == 1462


start_date = (2001, 01, 01, 0, 0)
end_date = (2003, 12, 31, 0, 0)
idxtimebeg, idxtimeend = ClimateTools.timeindex(d, start_date, end_date, T)
@test idxtimebeg == 367
@test idxtimeend == 1461

end
