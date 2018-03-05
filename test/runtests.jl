using ClimateTools
using AxisArrays
# using Lint
using Base.Test

# TODO tests for annualmean and annualsum
# @test isempty(lintpkg("ClimateTools"))

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
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = prcp1(C)
@test ind.data.data == Results

# annualsum
d = Date(2003,1,1):Date(2005,12,31)

data = vcat(ones(Float64, 365), zeros(Float64, 366), ones(Float64, 365))
Results = Array{Int64, 1}(3); Results[1] = 365.; Results[2] = 0.; Results[3] = 365.;
@test annualsum(data, d) == Results

data= vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 2}(3, 1); Results[1] = 365.; Results[2] = 0.; Results[3] = 365.;
@test annualsum(data, d) == Results

data = Array{Float64,3}(1096, 2, 2)
data[:,1,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,1,2] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,2] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(3, 2, 2); Results[1,1,1] = 365.; Results[2,1,1] = 0.; Results[3,1,1] = 365.; Results[1,2,1] = 365.; Results[2,2,1] = 0; Results[3,2,1] = 365.; Results[1,1,2] = 365.; Results[2,1,2] = 0.; Results[3,1,2] = 365.; Results[1,2,2] = 365.; Results[2,2,2] = 0.; Results[3,2,2] = 365.;
@test annualsum(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = annualsum(C)
@test ind.data.data == Results

# annualmean
d = Date(2003,1,1):Date(2005,12,31)

data = vcat(ones(Float64, 365), zeros(Float64, 366), ones(Float64, 365))
Results = Array{Int64, 1}(3); Results[1] = 1.; Results[2] = 0.; Results[3] = 1.;
@test annualmean(data, d) == Results

data= vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 2}(3, 1); Results[1] = 1.; Results[2] = 0.; Results[3] = 1.;
@test annualmean(data, d) == Results

data = Array{Float64,3}(1096, 2, 2)
data[:,1,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,1] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,1,2] = vcat(ones(365,1), zeros(366,1), ones(365))
data[:,2,2] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(3, 2, 2); Results[1,1,1] = 1.; Results[2,1,1] = 0.; Results[3,1,1] = 1.; Results[1,2,1] = 1.; Results[2,2,1] = 0; Results[3,2,1] = 1.; Results[1,1,2] = 1.; Results[2,1,2] = 0.; Results[3,1,2] = 1.; Results[1,2,2] = 1.; Results[2,2,2] = 0.; Results[3,2,2] = 1.;
@test annualmean(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = annualmean(C)
@test ind.data.data == Results

# Frostdays
data= vcat(-ones(365), -ones(366),zeros(365))
Results = Array{Int64, 1}(3); Results[1] = 365; Results[2] = 366; Results[3] = 0;
@test frostdays(data, d) == Results
@test icingdays(data, d) == Results


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
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = frostdays(C)
@test ind.data.data == Results
C = ClimateTools.ClimGrid(axisdata, variable = "tasmax")
ind = icingdays(C)
@test ind.data.data == Results

# Summer Days
d = Date(2003,1,1):Date(2007,12,31)
data = collect(1.:1826.)
@test summerdays(data, d) == [340, 366, 365, 365, 365]

data = hcat(collect(1.:1826.), collect(1.:1826.))
@test summerdays(data, d) == hcat([340, 366, 365, 365, 365],[340, 366, 365, 365, 365])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(1.:1826.); data[:, 1, 2] = collect(1.:1826.);data[:, 2, 1] = collect(1.:1826.);data[:, 2, 2] = collect(1.:1826.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [340, 366, 365, 365, 365]''; Results[:, 1, 2] = [340, 366, 365, 365, 365]''; Results[:, 2, 1] = [340, 366, 365, 365, 365]''; Results[:, 2, 2] = [340, 366, 365, 365, 365]'';
@test summerdays(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmax")
ind = summerdays(C)
@test ind.data.data == Results

# Tropical Nights
d = Date(2003,1,1):Date(2007,12,31)
data = collect(1.:1826.)
@test tropicalnights(data, d) == [345, 366, 365, 365, 365]

data = hcat(collect(1.:1826.), collect(1.:1826.))
@test tropicalnights(data, d) == hcat([345, 366, 365, 365, 365],[345, 366, 365, 365, 365])

data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(1.:1826.); data[:, 1, 2] = collect(1.:1826.);data[:, 2, 1] = collect(1.:1826.);data[:, 2, 2] = collect(1.:1826.);
Results = Array{Int64, 3}(5, 2, 2); Results[:, 1, 1] = [345, 366, 365, 365, 365]''; Results[:, 1, 2] = [345, 366, 365, 365, 365]''; Results[:, 2, 1] = [345, 366, 365, 365, 365]''; Results[:, 2, 2] = [345, 366, 365, 365, 365]'';
@test tropicalnights(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
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
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
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
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
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
Results = Array{Float64, 3}(5, 2, 2); Results[:, 1, 1] = [-436., -70., 295., 660., 1025.]''; Results[:, 1, 2] = [-436., -70., 295., 660., 1025.]''; Results[:, 2, 1] = [-436., -70., 295., 660., 1025.]''; Results[:, 2, 2] = [-436., -70., 295., 660., 1025.]'';
@test annualmax(data, d) == Results



# ClimGrid based tests
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = annualmax(C)
@test ind.data.data == Results

datetmp = Date(2003,01,01):Date(2008,1,1)
# idx = Dates.monthday(datetmp) .== (2,29)
idx = (Dates.month.(datetmp) .== 2) .& (Dates.day.(datetmp) .== 29)
datetmp = datetmp[.!idx] #creates an Array of dates
data = Array{Float64, 3}(1826, 2, 2)
data[:, 1, 1] = collect(-800.:1025.); data[:, 1, 2] = collect(-800.:1025.);data[:, 2, 1] = collect(-800.:1025.);data[:, 2, 2] = collect(-800.:1025.);
Results = Array{Float64, 3}(6, 2, 2); Results[:, 1, 1] = [-436., -71., 294., 659., 1024., 1025.]''; Results[:, 1, 2] = [-436., -71., 294., 659., 1024., 1025]''; Results[:, 2, 1] = [-436., -71., 294., 659., 1024., 1025.]''; Results[:, 2, 2] = [-436., -71., 294., 659., 1024., 1025.]'';
@test annualmax(data, datetmp) == Results

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
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = annualmin(C)
@test ind.data.data == Results

# Shapefile test
filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")

polyshp = read(filename,Shapefile.Handle)
x, y = shapefile_coords(polyshp.shapes[1])
P = [x y]
P = P'
lat = NetCDF.ncread(filenc, "lat")
lon = NetCDF.ncread(filenc, "lon")
msk = inpolyvec(lon, lat, P)
C = nc2julia(filenc, "tas", poly = P)
status, figh = mapclimgrid(C); @test status == true; PyPlot.close()
C = nc2julia(filenc, "tas")
status, figh = mapclimgrid(C, mask = msk); @test status == true; PyPlot.close()
@test length(x) == 6
@test length(y) == 6
@test isnan(x[1])
@test isnan(y[1])

filename = joinpath(dirname(@__FILE__), "data", "zoneAgricoleQc15km.shp")
polyshp = read(filename,Shapefile.Handle)
x, y = shapefile_coords(polyshp.shapes[2])
P = [x y]
P = P'
cities_lat = collect(45.5016889:0.1:46.8138783)
cities_lon = collect(-73.56725599999999:0.1:-71.2079809) + 360
@test sum(.!isnan.(inpolyvec(cities_lon, cities_lat, P))) == 336
@test sum(.!isnan.(inpolyvec(cities_lon - 360, cities_lat, P))) == 336

# Mapping test
filename = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
C = nc2julia(filename, "tas")
status, figh = mapclimgrid(C);@test status == true;PyPlot.close()
status, figh = mapclimgrid(C, region = "World");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Canada");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Quebec");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "NorthAmerica");@test status == true; PyPlot.close()
status, figh = mapclimgrid(annualmax(C), region = "Europe");@test status == true; PyPlot.close()

# precip
C = nc2julia(filename, "pr")
status, figh = mapclimgrid(C);@test status == true; PyPlot.close()
status, figh = mapclimgrid(prcp1(C), region = "World");@test status == true; PyPlot.close()
status, figh = mapclimgrid(prcp1(C), region = "Canada");@test status == true; PyPlot.close()
status, figh = mapclimgrid(prcp1(C), region = "Europe");@test status == true; PyPlot.close()
status, figh = mapclimgrid(prcp1(C), region = "NorthAmerica");@test status == true; PyPlot.close()

# ua wind
filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
polyshp = read(filename,Shapefile.Handle)
x, y = shapefile_coords(polyshp.shapes[1])
P = [x y]
P = P'
C = nc2julia(filenc, "ua")
status, figh = mapclimgrid(C, level = 3);@test status == true; PyPlot.close() # given level
status, figh = mapclimgrid(C);@test status == true; PyPlot.close() # feeding a 4D field
status, figh = mapclimgrid(C, poly = P);@test status == true; PyPlot.close() # feeding a 4D field with a polygon
status, figh = mapclimgrid(C, mask = msk);@test status == true; PyPlot.close() # feeding a 4D field with a mask
C = nc2julia(filenc, "ua", poly = P)
status, figh = mapclimgrid(C);@test status == true; PyPlot.close() # feeding a 4D field

# test that nc2julia return a ClimGrid type
filename = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
C = nc2julia(filename, "tas")
@test nc2julia(filename, "tas", data_units = "Celsius")[2] == "Celsius"
@test nc2julia(filename, "pr", data_units = "mm")[2] == "mm"
@test typeof(nc2julia(filename, "tas")) == ClimateTools.ClimGrid{AxisArrays.AxisArray{Float32,3,Array{Float32,3},Tuple{AxisArrays.Axis{:time,Array{Date,1}},AxisArrays.Axis{:lon,Array{Float32,1}},AxisArrays.Axis{:lat,Array{Float32,1}}}}}

@test typeof(buildtimevec(filename)) == Array{Date, 1}

# Time units
units = NetCDF.ncgetatt(filename, "time", "units") # get starting date
m = match(r"(\d+)[-.\/](\d+)[-.\/](\d+)", units, 1) # match a date from string
daysfrom = m.match # get only the date ()"yyyy-mm-dd" format)
initDate = Date(daysfrom, "yyyy-mm-dd")
timeRaw = floor.(NetCDF.ncread(filename, "time"))
@test sumleapyear(initDate::Date, timeRaw) == 485

# INTERFACE
B = vcat(C, C)
@test size(B.data) == (2, 256, 128) # vcat does not look at dimensions
B = merge(C, C)
@test size(B.data) == (1, 256, 128) # C being similar, they should not add up, as opposed to vcat
# Operators +, -, *, /
B = C + C
@test B[1].data[1, 1, 1] == 431.787f0
B = C - C
@test B[1].data[1, 1, 1] == 0.0
B = C / 2
@test B[1].data[1, 1, 1] == 107.94675f0
B = C / 2.2
@test B[1].data[1, 1, 1] == 98.13340620561078
B = C * 2
@test B[1].data[1, 1, 1] == 431.787f0
B = C * 2.2
@test B[1].data[1, 1, 1] == 474.9656860351563

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
@test_throws ErrorException C[12]
@test annualmax(C)[10] == "tas"
@test size(C) == (11, )
@test size(C, 1) == 11
@test length(C) == 11
@test endof(C) == 11
@test C[end] == "noleap"
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
C = nc2julia(filenc, "tas")
Csub = temporalsubset(C, Date(2000, 05, 15), Date(2000, 05, 15))
@test Csub[1][1, 1, 1] == 215.8935f0
@test Csub[1][Axis{:time}][1] == Date(2000, 05, 15)

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

# Test Interpolation
filename = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
C = nc2julia(filename, "tas")
# Get lat lon vector
lat = C[1][Axis{:lat}][:]
lon = C[1][Axis{:lon}][:]
# Shift longitude by 1
lon += 1
axisdata = AxisArray(C[1].data, Axis{:time}(C[1][Axis{:time}][:]), Axis{:lon}(lon), Axis{:lat}(lat))
C2 = ClimGrid(axisdata, variable = "tas")
@test interp_climgrid(C, C2)[1].data[1, 1, 1] == 215.83078002929688
@test interp_climgrid(C, lon, lat)[1].data[1, 1, 1] == 215.83078002929688

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
