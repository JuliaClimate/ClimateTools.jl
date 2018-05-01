# test that nc2julia return a ClimGrid type
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
C = nc2julia(filenc, "tas")
@test nc2julia(filenc, "tas", data_units = "Celsius")[2] == "Celsius"
@test nc2julia(filenc, "pr", data_units = "mm")[2] == "mm"
@test typeof(nc2julia(filenc, "tas")) == ClimateTools.ClimGrid{AxisArrays.AxisArray{Float32,3,Array{Float32,3},Tuple{AxisArrays.Axis{:time,Array{Date,1}},AxisArrays.Axis{:lon,Array{Float32,1}},AxisArrays.Axis{:lat,Array{Float32,1}}}}}

@test typeof(buildtimevec(filenc)) == Array{Date, 1}

# Time units
units = NetCDF.ncgetatt(filenc, "time", "units") # get starting date
m = match(r"(\d+)[-.\/](\d+)[-.\/](\d+)", units, 1) # match a date from string
daysfrom = m.match # get only the date ()"yyyy-mm-dd" format)
initDate = Date(daysfrom, "yyyy-mm-dd")
timeRaw = floor.(NetCDF.ncread(filenc, "time"))
@test sumleapyear(initDate::Date, timeRaw) == 485

# INTERFACE
B = vcat(C, C)
@test size(B.data) == (2, 256, 128) # vcat does not look at dimensions
B = merge(C, C)
@test size(B.data) == (1, 256, 128) # C being similar, they should not add up, as opposed to vcat
# Operators +, -, *, /
B = C + C
@test B[1].data[1, 1, 1] == 438.4457f0
B = C - C
@test B[1].data[1, 1, 1] == 0.0
B = C / 2
@test B[1].data[1, 1, 1] == 109.61143f0
B = C / 2.2
@test B[1].data[1, 1, 1] == 99.6467520973899
B = C * 2
@test B[1].data[1, 1, 1] == 438.4457f0
B = C * 2.2
@test B[1].data[1, 1, 1] == 482.2902801513672

# @test typeof(show(C)) == Dict{Any, Any}
@test typeof(C[1].data) == Array{Float64,3} || typeof(C[1].data) == Array{Float32,3}
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
@test size(C) == (21, )
@test size(C, 1) == 21
@test length(C) == 21
@test endof(C) == 21
@test_throws ErrorException C[end]
@test ndims(C) == 1


# Spatial subset
filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
polyshp = read(filename,Shapefile.Handle)
x, y = shapefile_coords(polyshp.shapes[1])
P = [x y]
P = P'
C = nc2julia(filenc, "tas")
Csub = spatialsubset(C, P)
@test size(Csub[1]) == (1, 23, 12)
@test Csub[1][1, 1, 1] == 294.6609f0
Csub = spatialsubset(C, P')
@test size(Csub[1]) == (1, 23, 12)
@test Csub[1][1, 1, 1] == 294.6609f0
C = nc2julia(filenc, "ua")
Csub = spatialsubset(C, P)
@test size(Csub[1]) == (1, 23, 12, 17)
@test Csub[1][1, 12, 1, 1] == 6.658482f0
@test isnan(Csub[1][1, 1, 1, 1])

poly= [[NaN 10 -10 -10 10 10];[NaN -10 -20 10 10 -10]] # meridian test
C = nc2julia(filenc, "tas", poly=poly)

# Spatial subset
C = nc2julia(filenc, "tas")
Csub = temporalsubset(C, Date(2000, 05, 15), Date(2000, 05, 15))
@test Csub[1][1, 1, 1] == 219.22285f0
@test Csub[1][Axis{:time}][1] == Date(2000, 05, 15)
B = nc2julia(filenc, "tas", start_date=Date(2000, 05, 15), end_date=Date(2000, 05, 15))
@test B[1] == C[1]

# Time resolution
timevec = [1, 2, 3]
@test timeresolution(timevec) == "24h"
timevec = [1.0, 1.5, 2.0]
@test timeresolution(timevec) == "12h"
timevec = [1.25, 1.5, 1.75]
@test timeresolution(timevec) == "6h"
timevec = [1.125, 1.25, 1.375]
@test timeresolution(timevec) == "3h"
timevec = NetCDF.ncread(filenc, "time")
@test timeresolution(timevec) == "N/A"



# MESHGRID
YV = [1, 2, 3]
XV = [1, 2, 3]
@test meshgrid(XV, YV) == ([1 2 3; 1 2 3; 1 2 3], [1 1 1; 2 2 2; 3 3 3])

## INPOLY
@test leftorright(0.5,0.5, 1,0,1,1) == -1
@test leftorright(1.5,.5, 1,0,1,1) == 1
@test leftorright(1,0.5, 1,0,1,1) == 0

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
if VERSION >= v"0.5-"
    eval(:(@test_broken inpoly(p1, poly) )) # should be true
end
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

# Test Interpolation
filename = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
C = nc2julia(filename, "tas")
# Get lat lon vector
lat = Float32.(C[1][Axis{:lat}][:])
lon = Float32.(C[1][Axis{:lon}][:])
latgrid = Float32.(C.latgrid)
longrid = Float32.(C.longrid)
# Shift longitude by 1
lon += Float32(1.0)
longrid += Float32(1.0)
axisdata = AxisArray(C[1].data, Axis{:time}(C[1][Axis{:time}][:]), Axis{:lon}(lon), Axis{:lat}(lat))
C2 = ClimGrid(axisdata, variable = "tas", longrid=longrid, latgrid=latgrid, msk=C.msk)
@test interp_climgrid(C, C2)[1].data[1, 1, 1] == 219.2400638156467
@test interp_climgrid(C, lon, lat)[1].data[1, 1, 1] == 219.2400638156467

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
data = randn(3, 3, 2)
mask = [[NaN; 1; 1] [1.; NaN;1.]]
@test isnan(applymask(data, mask)[1, 1, 1]) && isnan(applymask(data, mask)[1, 2, 2]) && isnan(applymask(data, mask)[2, 1, 1]) && isnan(applymask(data, mask)[2, 2, 2]) && isnan(applymask(data, mask)[3, 1, 1]) && isnan(applymask(data, mask)[3, 2, 2])

for i = 1:size(data, 1)
    @test applymask(data, mask)[i, 2, 1] == data[i, 2, 1]
    @test applymask(data, mask)[i, 1, 2] == data[i, 1, 2]
    @test applymask(data, mask)[i, 3, 1] == data[i, 3, 1]
    @test applymask(data, mask)[i, 3, 2] == data[i, 3, 2]
end

# 4-D data
data = randn(3, 3, 2, 1)
mask = [[NaN; 1; 1] [1.; NaN;1.]]
@test isnan(applymask(data, mask)[1, 1, 1, 1]) && isnan(applymask(data, mask)[1, 2, 2, 1]) && isnan(applymask(data, mask)[2, 1, 1, 1]) && isnan(applymask(data, mask)[2, 2, 2, 1]) && isnan(applymask(data, mask)[3, 1, 1, 1]) && isnan(applymask(data, mask)[3, 2, 2, 1])

for i = 1:size(data, 1)
    @test applymask(data, mask)[i, 2, 1, 1] == data[i, 2, 1, 1]
    @test applymask(data, mask)[i, 1, 2, 1] == data[i, 1, 2, 1]
    @test applymask(data, mask)[i, 3, 1, 1] == data[i, 3, 1, 1]
    @test applymask(data, mask)[i, 3, 2, 1] == data[i, 3, 2, 1]
end

# Test sumleapyear with StepRange{Date,Base.Dates.Day} type
d = Date(2003,1,1):Date(2008,12,31)
@test sumleapyear(d) == 2

# Test timeresolution and pr_timefactor
filename = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
timevec = NetCDF.ncread(filename, "time")
@test pr_timefactor(timeresolution(timevec)) == 1.
@test pr_timefactor("24h") == 86400.
@test pr_timefactor("12h") == 43200.
@test pr_timefactor("6h") == 21600.
@test pr_timefactor("3h") == 10800.
