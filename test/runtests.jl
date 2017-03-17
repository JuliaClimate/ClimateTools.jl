using ClimateTools
using AxisArrays
using Lint
using Base.Test

@test isempty(lintpkg("ClimateTools"))

# Prcp1
d = Date(2003,1,1):Date(2005,12,31)

data = vcat(ones(Float64, 365), zeros(Float64, 366), ones(Float64, 365))
Results = Array{Int64, 1}(3); Results[1] = 365; Results[2] = 0; Results[3] = 365;
@test prcp1(data, d) == Results

data= vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 2}(3, 1); Results[1] = 365; Results[2] = 0; Results[3] = 365;
@test prcp1(data, d) == Results

data = Array{Float64,3}(1096, 2, 2)
data[:,1,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,1,2] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,2] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(3, 2, 2); Results[1,1,1] = 365; Results[2,1,1] = 0; Results[3,1,1] = 365; Results[1,2,1] = 365; Results[2,2,1] = 0; Results[3,2,1] = 365; Results[1,1,2] = 365; Results[2,1,2] = 0; Results[3,1,2] = 365; Results[1,2,2] = 365; Results[2,2,2] = 0; Results[3,2,2] = 365;
@test prcp1(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, var = "pr")
ind = prcp1(C)
@test ind.data.data == Results

# Frostdays
data= vcat(-ones(365,1), -ones(366),zeros(365,1))
Results = Array{Int64, 2}(3, 1); Results[1] = 365; Results[2] = 366; Results[3] = 0;
@test frostdays(data, d) == Results
@test icingdays(data, d) == Results

data = Array{Float64,3}(1096, 2, 2)
data[:,1,1] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[:,2,1] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[:,1,2] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[:,2,2] = vcat(-ones(365,1), -ones(366,1), zeros(365))
Results = Array{Int64, 3}(3, 2, 2); Results[1,1,1] = 365; Results[2,1,1] = 366; Results[3,1,1] = 0; Results[1,2,1] = 365; Results[2,2,1] = 366; Results[3,2,1] = 0; Results[1,1,2] = 365; Results[2,1,2] = 366; Results[3,1,2] = 0; Results[1,2,2] = 365; Results[2,2,2] = 366; Results[3,2,2] = 0;
@test frostdays(data, d) == Results
@test icingdays(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, var = "tasmin")
ind = frostdays(C)
@test ind.data.data == Results
C = ClimateTools.ClimGrid(axisdata, var = "tasmax")
ind = icingdays(C)
@test ind.data.data == Results

# Summer Days
d = Date(2003,1,1):Date(2007,12,31)
data = collect(1.:1826.)
@test summerdays(data, d) == [341, 366, 365, 365, 365]

data = hcat(collect(1.:1826.), collect(1.:1826.))
@test summerdays(data, d) == hcat([341, 366, 365, 365, 365],[341, 366, 365, 365, 365])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(1.:1826.); data[:, 1, 2] = collect(1.:1826.);data[:, 2, 1] = collect(1.:1826.);data[:, 2, 2] = collect(1.:1826.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [341, 366, 365, 365, 365]''; Results[:, 1, 2] = [341, 366, 365, 365, 365]''; Results[:, 2, 1] = [341, 366, 365, 365, 365]''; Results[:, 2, 2] = [341, 366, 365, 365, 365]'';
@test summerdays(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, var = "tasmax")
ind = summerdays(C)
@test ind.data.data == Results

# Tropical Nights
d = Date(2003,1,1):Date(2007,12,31)
data = collect(1.:1826.)
@test tropicalnights(data, d) == [346, 366, 365, 365, 365]

data = hcat(collect(1.:1826.), collect(1.:1826.))
@test tropicalnights(data, d) == hcat([346, 366, 365, 365, 365],[346, 366, 365, 365, 365])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(1.:1826.); data[:, 1, 2] = collect(1.:1826.);data[:, 2, 1] = collect(1.:1826.);data[:, 2, 2] = collect(1.:1826.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [346, 366, 365, 365, 365]''; Results[:, 1, 2] = [346, 366, 365, 365, 365]''; Results[:, 2, 1] = [346, 366, 365, 365, 365]''; Results[:, 2, 2] = [346, 366, 365, 365, 365]'';
@test tropicalnights(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, var = "tasmin")
ind = tropicalnights(C)
@test ind.data.data == Results

# Custom thresholds
d = Date(2003,1,1):Date(2007,12,31)
data = collect(1.:1826.)
@test customthresover(data, d, 20) == [345, 366, 365, 365, 365]

data = hcat(collect(1.:1826.), collect(1.:1826.))
@test customthresover(data, d, 20) == hcat([345, 366, 365, 365, 365],[345, 366, 365, 365, 365])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(1.:1826.); data[:, 1, 2] = collect(1.:1826.);data[:, 2, 1] = collect(1.:1826.);data[:, 2, 2] = collect(1.:1826.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [345, 366, 365, 365, 365]''; Results[:, 1, 2] = [345, 366, 365, 365, 365]''; Results[:, 2, 1] = [345, 366, 365, 365, 365]''; Results[:, 2, 2] = [345, 366, 365, 365, 365]'';
@test customthresover(data, d, 20) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, var = "tasmin")
ind = customthresover(C, 20)
@test ind.data.data == Results

d = Date(2003,1,1):Date(2007,12,31)
data = collect(-800.:1025.)
@test customthresunder(data, d, 200) == [365, 366, 269, 0, 0]

data = hcat(collect(-800.:1025.), collect(-800.:1025.))
@test customthresunder(data, d, 200) == hcat([365, 366, 269, 0, 0],[365, 366, 269, 0, 0])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(-800.:1025.); data[:, 1, 2] = collect(-800.:1025.);data[:, 2, 1] = collect(-800.:1025.);data[:, 2, 2] = collect(-800.:1025.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [365, 366, 269, 0, 0]''; Results[:, 1, 2] = [365, 366, 269, 0, 0]''; Results[:, 2, 1] = [365, 366, 269, 0, 0]''; Results[:, 2, 2] = [365, 366, 269, 0, 0]'';
@test customthresunder(data, d, 200) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, var = "tasmin")
ind = customthresunder(C, 200)
@test ind.data.data == Results

# ANNUAL MAXIMUM
d = Date(2003,1,1):Date(2007,12,31)
data = collect(-800.:1025.)
@test annualmax(data, d) == [-436., -70., 295., 660., 1025.]

data = hcat(collect(-800.:1025.), collect(-800.:1025.))
@test annualmax(data, d) == hcat([-436., -70., 295., 660., 1025.],[-436., -70., 295., 660., 1025.])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(-800.:1025.); data[:, 1, 2] = collect(-800.:1025.);data[:, 2, 1] = collect(-800.:1025.);data[:, 2, 2] = collect(-800.:1025.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [-436., -70., 295., 660., 1025.]''; Results[:, 1, 2] = [-436., -70., 295., 660., 1025.]''; Results[:, 2, 1] = [-436., -70., 295., 660., 1025.]''; Results[:, 2, 2] = [-436., -70., 295., 660., 1025.]'';
@test annualmax(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, var = "tasmin")
ind = annualmax(C)
@test ind.data.data == Results

# ANNUAL MINIMUM
d = Date(2003,1,1):Date(2007,12,31)
data = collect(-800.:1025.)
@test annualmin(data, d) == [-800., -435., -69., 296., 661.]

data = hcat(collect(-800.:1025.), collect(-800.:1025.))
@test annualmin(data, d) == hcat([-800., -435., -69., 296., 661.],[-800., -435., -69., 296., 661.])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(-800.:1025.); data[:, 1, 2] = collect(-800.:1025.);data[:, 2, 1] = collect(-800.:1025.);data[:, 2, 2] = collect(-800.:1025.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [-800., -435., -69., 296., 661.]''; Results[:, 1, 2] = [-800., -435., -69., 296., 661.]''; Results[:, 2, 1] = [-800., -435., -69., 296., 661.]''; Results[:, 2, 2] = [-800., -435., -69., 296., 661.]'';
@test annualmin(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, var = "tasmin")
ind = annualmin(C)
@test ind.data.data == Results

# # Mapping test
d = Date(2003,1,1):Date(2005,12,31)
data = randn(1096, 200, 81)
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(101:300), Axis{:lat}(-20:60))
C = ClimateTools.ClimGrid(axisdata, var = "pr")
@test mapclimgrid(C)
# @test mapclimgrid(prcp1(C)) == 1
# @test mapclimgrid(prcp1(C), region = "World") == 1
# @test mapclimgrid(prcp1(C), region = "Canada") == 1
# C = ClimateTools.ClimGrid(axisdata, var = "tasmax")
# @test mapclimgrid(C) == 1
# @test mapclimgrid(annualmax(C)) == 1
# @test mapclimgrid(annualmax(C), region = "World") == 1
# @test mapclimgrid(annualmax(C), region = "Canada") == 1

# NetCDF Extraction test
filename = joinpath(Pkg.dir("ClimateTools"), "test", "data", "sresa1b_ncar_ccsm3-example.nc")
# nc2julia(filename)
C = nc2julia(filename, "tas", poly = [0. 0.])
@test typeof(nc2julia(filename, "tas", poly = [0. 0.])) == ClimateTools.ClimGrid
@test typeof(buildtimevec(filename)) == Array{Date, 1}
# Time units
units = NetCDF.ncgetatt(filename, "time", "units") # get starting date
m = match(r"(\d+)[-.\/](\d+)[-.\/](\d+)", units, 1) # match a date from string
daysfrom = m.match # get only the date ()"yyyy-mm-dd" format)
initDate = Date(daysfrom, "yyyy-mm-dd")
timeRaw = floor(NetCDF.ncread(filename, "time"))
@test sumleapyear(initDate::Date, timeRaw) == 485

# INTERFACE
@test typeof(vcat(C, C)) == ClimateTools.ClimGrid
# @test typeof(show(C)) == Dict{Any, Any}
@test typeof(C[1].data) == Array{Float64,3} || Array{Float32,3}
@test C[2] == "Celsius"
@test C[3] == ""
@test C[4] == "720 ppm stabilization experiment (SRESA1B)"
@test C[5] == ""
@test C[6] == "degrees_east"
@test C[7] == "degrees_north"
@test C[8] == joinpath(Pkg.dir("ClimateTools"), "test", "data", "sresa1b_ncar_ccsm3-example.nc")
@test C[9] == "tas"
@test annualmax(C)[9] == "annualmax"
@test C[10] == "tas"
@test annualmax(C)[10] == "tas"
@test size(C) == (10, )
@test size(C, 1) == 10
@test length(C) == 10

# MESHGRID
YV = [1 2 3]'
XV = [1 2 3]'
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
