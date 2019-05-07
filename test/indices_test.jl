@testset "Indices" begin

# Prcp1
d = collect(DateTime(2003,1,1):Day(1):DateTime(2005,12,31))

data = Array{Float64,3}(undef, 2, 2,1096)
data[1,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[1,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(undef, 2, 2, 3);
Results[1,1,1] = 365;
Results[1,1, 2] = 0;
Results[1,1,3] = 365;
Results[1,2,1] = 365;
Results[2,1,1] = 365;
Results[2,1,2] = 0;
Results[2,1,3] = 365;
Results[1,2,1] = 365;
Results[1,2,2] = 0;
Results[1,2,3] = 365;
Results[2,2,1] = 365;
Results[2,2,2] = 0;
Results[2,2,3] = 365;

# # Add units
# data = [data][1]u"mm"

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = prcp1(C)
@test ind.data.data == Results

# annualsum
data = Array{Float64,3}(undef, 2, 2, 1096)
data[1,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[1,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(undef, 2, 2, 3);
Results[1,1,1] = 365.;
Results[1,1,2] = 0.;
Results[1,1,3] = 365.;
Results[2,1,1] = 365.;
Results[2,1,2] = 0;
Results[2,1,3] = 365.;
Results[1,2,1] = 365.;
Results[1,2,2] = 0.;
Results[1,2,3] = 365.;
Results[2,2,1] = 365.;
Results[2,2,2] = 0.;
Results[2,2,3] = 365.;

# Add units
# data = [data][1]u"mm"
# Results = [Results][1]u"mm"

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = annualsum(C)
@test ind.data.data == Results

# annualmean
data = Array{Float64,3}(undef, 2, 2, 1096)
data[1,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,1,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[1,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
data[2,2,:] = vcat(ones(365,1), zeros(366,1), ones(365))
Results = Array{Int64, 3}(undef, 2, 2, 3);
Results[1,1,1] = 1.;
Results[1,1,2] = 0.;
Results[1,1,3] = 1.;
Results[2,1,1] = 1.;
Results[2,1,2] = 0;
Results[2,1,3] = 1.;
Results[1,2,1] = 1.;
Results[1,2,2] = 0.;
Results[1,2,3] = 1.;
Results[2,2,1] = 1.;
Results[2,2,2] = 0.;
Results[2,2,3] = 1.;

# Add units
# data = [data][1]u"mm"
# Results = [Results][1]u"mm"

# @test annualmean(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "pr")
ind = annualmean(C)
@test ind.data.data == Results

# Frostdays

data = Array{Float64,3}(undef, 2, 2, 1096)
data[1,1,:] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[2,1,:] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[1,2,:] = vcat(-ones(365,1), -ones(366,1), zeros(365))
data[2,2,:] = vcat(-ones(365,1), -ones(366,1), zeros(365))
Results = Array{Int64, 3}(undef, 2, 2, 3);
Results[1,1,1] = 365;
Results[1,1,2] = 366;
Results[1,1,3] = 0;
Results[2,1,1] = 365;
Results[2,1,2] = 366;
Results[2,1,3] = 0;
Results[1,2,1] = 365;
Results[1,2,2] = 366;
Results[1,2,3] = 0;
Results[2,2,1] = 365;
Results[2,2,2] = 366;
Results[2,2,3] = 0;

# Add units
# data = [data][1]u"°C"

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = frostdays(C)
@test ind.data.data == Results
C = ClimateTools.ClimGrid(axisdata, variable = "tasmax")
ind = icingdays(C)
@test ind.data.data == Results

# Summer Days
d = collect(DateTime(2003,1,1):Day(1):DateTime(2007,12,31))

data = Array{Float64, 3}(undef, 2, 2, 1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:] = collect(1.0:1826.0); data[2,1,:]=collect(1.0:1826.0); data[2,2,:] = collect(1.0:1826.0);
Results = Array{Int64, 3}(undef, 2, 2, 5);
Results[1, 1, :]=[340, 366, 365, 365, 365]'';
Results[1, 2, :] = [340, 366, 365, 365, 365]'';
Results[2, 1, :] = [340, 366, 365, 365, 365]'';
Results[2, 2, :] = [340, 366, 365, 365, 365]'';

# Add units
# data = [data][1]u"°C"

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmax")
ind = summerdays(C)
@test ind.data.data == Results

data = Array{Float64, 3}(undef, 2, 2, 1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:]=collect(1.0:1826.0); data[2,1,:]=collect(1.0:1826.0); data[2,2,:]=collect(1.0:1826.0);
Results = Array{Int64, 3}(undef, 2,2,5);
Results[1,1,:]=[345, 366, 365, 365, 365]'';
Results[1,2,:]=[345, 366, 365, 365, 365]'';
Results[2,1,:]=[345, 366, 365, 365, 365]'';
Results[2,2,:]=[345, 366, 365, 365, 365]'';

# Add units
# data = [data][1]u"°C"

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = tropicalnights(C)
@test ind.data.data == Results

# Custom thresholds

data = Array{Float64, 3}(undef, 2,2,1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:]=collect(1.0:1826.0); data[2,1,:]=collect(1.0:1826.0); data[2,2,:] = collect(1.0:1826.0);
Results = Array{Int64, 3}(undef, 2,2,5);
Results[1,1,:] = [345, 366, 365, 365, 365]'';
Results[1,2,:] = [345, 366, 365, 365, 365]'';
Results[2,1,:] = [345, 366, 365, 365, 365]'';
Results[2,2,:] = [345, 366, 365, 365, 365]'';

# Add units
# data = [data][1]u"°C"

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = customthresover(C, 20)
@test ind.data.data == Results

data = Array{Float64, 3}(undef, 2,2,1826)
data[1,1,:] = collect(-800.0:1025.0); data[1,2,:] = collect(-800.0:1025.0); data[2,1,:]=collect(-800.0:1025.0); data[2,2,:] = collect(-800.0:1025.0);
Results = Array{Int64, 3}(undef, 2,2,5);
Results[1,1,:]=[365, 366, 269, 0, 0]'';
Results[1,2,:]=[365, 366, 269, 0, 0]'';
Results[2,1,:]=[365, 366, 269, 0, 0]'';
Results[2,2,:]=[365, 366, 269, 0, 0]'';

# # Add units
# data = [data][1]u"°C"

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = customthresunder(C, 200)
@test ind.data.data == Results

# ANNUAL MAXIMUM
data = Array{Float64, 3}(undef, 2,2,1826)
data[1,1,:] = collect(-800.0:1025.0); data[1,2,:] = collect(-800.0:1025.0); data[2,1,:] = collect(-800.0:1025.0); data[2,2,:] = collect(-800.0:1025.0);
Results = Array{Float64, 3}(undef, 2,2,5);
Results[1,1,:]=[-436.0, -70.0, 295.0, 660.0, 1025.0]'';
Results[1,2,:]=[-436.0, -70.0, 295.0, 660.0, 1025.0]'';
Results[2,1,:]=[-436.0, -70.0, 295.0, 660.0, 1025.0]'';
Results[2,2,:]=[-436.0, -70.0, 295.0, 660.0, 1025.0]'';

# # Add units
# data = [data][1]u"K"
# Results = [Results][1]u"K"

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = annualmax(C)
@test ind.data.data == Results

# ANNUAL MINIMUM
data = Array{Float64, 3}(undef, 2,2,1826)
data[1,1,:] = collect(-800.:1025.); data[1,2,:]=collect(-800.:1025.); data[2,1,:]=collect(-800.:1025.); data[2,2,:]=collect(-800.:1025.);
Results = Array{Int64, 3}(undef, 2,2,5);
Results[1,1,:]=[-800.0,-435.0,-69.0,296.0,661.0]'';
Results[1,2,:]=[-800.0,-435.0,-69.0,296.0,661.0]'';
Results[2,1,:]=[-800.0,-435.0,-69.0,296.0,661.0]'';
Results[2,2,:]=[-800.0,-435.0,-69.0,296.0,661.0]'';
 #
 # # Add units
 # data = [data][1]u"K"
 # Results = [Results][1]u"K"

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = annualmin(C)
@test ind.data.data == Results


# Period mean
data = Array{Float64, 3}(undef, 2, 2, 1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:] = collect(-10.0:1815.0); data[2,1,:] = collect(0.0:1825.0); data[2,2,:] = collect(-1725.0:100);
Results = Array{Float64, 2}(undef, 2, 2);
Results[1, 1]= 746.0;
Results[1, 2] = 735.0;
Results[2, 1] = 745.0;
Results[2, 2] = -980.0;

# # Add units
# data = [data][1]u"K"
# Results = [Results][1]u"K"

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata)
ind = periodmean(C, start_date=(2004, 06, 01), end_date=(2005, 08, 31)) #mean between June 1st 2004 and August 31st 2005
@test ind.data.data == Results


Results = Array{Float64, 2}(undef, 2, 2);
Results[1, 1]= 913.5;
Results[1, 2] = 902.5;
Results[2, 1] = 912.5;
Results[2, 2] = -812.5;

# Add units
# Results = [Results][1]u"K"

ind = periodmean(C)
@test ind.data.data == Results

# Daymean
raw = collect(0:0.25:1459.75) # 6hr data
d = NCDatasets.timedecode(raw, "days since 2000-01-01", "noleap")

data = Array{Float64}(undef, 2, 2, 5840)
data[1, 1, :] = collect(1:5840)
data[1, 2, :] = collect(1:5840)
data[2, 1, :] = collect(1:5840)
data[2, 2, :] = collect(1:5840)

axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, frequency="6h")
D = daymean(C)
@test D[1][1,1,1] == 2.5
@test D[1][1,1,end] == 5838.5

D = daysum(C)
@test D[1][1,1,1] == 10.0
@test D[1][1,1,end] == 23354.0

@test get_timevec(C)[1] == DateTimeNoLeap(2000, 01, 01)


# Test vaporpressure(specific_humidity::ClimGrid, surface_pressure::ClimGrid)
d = DateTime(2003,1,1):Day(1):DateTime(2003,1,3)
# Dummy data
data_huss = Array{Float64,3}(undef, 2,2,3)
data_huss[1,1,:] .= 0
data_huss[1,2,:] .= 0.1
data_huss[2,1,:] .= 0.5
data_huss[2,2,:] .= 1
data_ps = Array{Float64,3}(undef, 2,2,3)
data_ps[:,:,1] .= 10.0
data_ps[:,:,2] .= 100.0
data_ps[:,:,3] .= 1000.0
# Expected values
Results = Array{Float64,3}(undef, 2,2,3)
Results[1,1,1] = 0.0;
Results[1,1,2] = 0.0;
Results[1,1,3] = 0.0;
Results[2,1,1] = data_huss[2,1,1]*data_ps[2,1,1]/(0.622+data_huss[2,1,1]);
Results[2,1,2] = data_huss[2,1,2]*data_ps[2,1,2]/(0.622+data_huss[2,1,2]);
Results[2,1,3] = data_huss[2,1,3]*data_ps[2,1,3]/(0.622+data_huss[2,1,3]);
Results[1,2,1] = data_huss[1,2,1]*data_ps[1,2,1]/(0.622+data_huss[1,2,1]);
Results[1,2,2] = data_huss[1,2,2]*data_ps[1,2,2]/(0.622+data_huss[1,2,2]);
Results[1,2,3] = data_huss[1,2,3]*data_ps[1,2,3]/(0.622+data_huss[1,2,3]);
Results[2,2,1] = data_huss[2,2,1]*data_ps[2,2,1]/(0.622+data_huss[2,2,1]);
Results[2,2,2] = data_huss[2,2,2]*data_ps[2,2,2]/(0.622+data_huss[2,2,2]);
Results[2,2,3] = data_huss[2,2,3]*data_ps[2,2,3]/(0.622+data_huss[2,2,3]);
# Using the function on dummy data
axisdata_huss = AxisArray(data_huss, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_ps = AxisArray(data_ps, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C_huss = ClimateTools.ClimGrid(axisdata_huss, variable = "huss")
C_ps = ClimateTools.ClimGrid(axisdata_ps, dataunits="Pa", variable = "ps")
vp = vaporpressure(C_huss, C_ps)
# Run the test
@test vp.data.data == Results


# Test  vaporpressure(specific_humidity::ClimGrid, sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)
# Dummy data
data_huss = Array{Float64,3}(undef, 2,2,3)
data_huss[1,1,:] .= 0
data_huss[1,2,:] .= 0.1
data_huss[2,1,:] .= 0.5
data_huss[2,2,:] .= 1

data_psl = Array{Float64,3}(undef, 2,2,3)
data_psl[:,:,1] .= 10.0
data_psl[:,:,2] .= 100.0
data_psl[:,:,3] .= 1000.0
# data_psl = [data_psl][1]u"Pa"

data_orog = Array{Float64,2}(undef, 2,2)
data_orog[:,1] .= 0.0
data_orog[:,2] .= 10.0
# data_orog = [data_orog][1]u"m"

data_tas = Array{Float64,3}(undef, 2,2,3)
data_tas[1,:,:] .= 100.0
data_tas[2,:,:] .= 300.0
# data_tas = [data_tas][1]u"K"

# Expected results
Results = Array{Float64,3}(undef, 2,2,3)
Results[1,1,:] .= round(0.0, digits=10)
Results[1,2,1] = round(1.380315267090330, digits=10)
Results[1,2,2] = round(13.80315267090330, digits=10)
Results[1,2,3] = round(138.0315267090330, digits=10)
Results[2,1,1] = round(4.456327985739750, digits=10)
Results[2,1,2] = round(44.56327985739750, digits=10)
Results[2,1,3] = round(445.6327985739750, digits=10)
Results[2,2,1] = round(6.158207427095860, digits=10)
Results[2,2,2] = round(61.58207427095860, digits=10)
Results[2,2,3] = round(615.8207427095860, digits=10)
# Results = [Results][1]u"Pa"

# Create the Climgrids
axisdata_huss = AxisArray(data_huss, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_psl = AxisArray(data_psl, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_tas = AxisArray(data_tas, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_psl = AxisArray(data_psl, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_orog = AxisArray(data_orog, Axis{:lon}(1:2), Axis{:lat}(1:2))
C_huss = ClimateTools.ClimGrid(axisdata_huss, variable = "huss")
C_psl = ClimateTools.ClimGrid(axisdata_psl, dataunits = "Pa", variable = "psl")
C_tas = ClimateTools.ClimGrid(axisdata_tas, dataunits = "K", variable = "tas")
C_orog = ClimateTools.ClimGrid(axisdata_orog, dataunits = "m", variable = "orog")
# Using the function on dummy data
vp = vaporpressure(C_huss, C_psl, C_orog, C_tas)
# Run the test
@test round.(vp.data.data, digits=10) == Results

# Test for tas in Celsius
data_tas = Array{Float64,3}(undef, 2,2,3)
data_tas[1,:,:] .= 100-273.15
data_tas[2,:,:] .= 300-273.15
# data_tas = [data_tas][1]u"°C"

axisdata_tas = AxisArray(data_tas, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C_tas = ClimateTools.ClimGrid(axisdata_tas, dataunits = "Celsius", variable = "tas")
# Using the function on dummy data
vp = vaporpressure(C_huss, C_psl, C_orog, C_tas)
# Run the test
@test round.(vp.data.data, digits=10) == Results

# Test wbgt(diurnal_temperature::ClimGrid, vapor_pressure::ClimGrid)
d = DateTime(2003,1,1):Day(1):DateTime(2003,1,3)
# Dummy data
data_tdiu = Array{Float64,3}(undef, 2,2,3)
data_tdiu[1,1,:] .= 250.0
data_tdiu[1,2,:] .= 275.0
data_tdiu[2,1,:] .= 300.0
data_tdiu[2,2,:] .= 325.0
# data_tdiu = [data_tdiu][1]u"K"

data_vp = Array{Float64,3}(undef, 2,2,3)
data_vp[:,:,1] .= 10.0
data_vp[:,:,2] .= 100.0
data_vp[:,:,3] .= 1000.0
# data_vp = [data_vp][1]u"Pa"

# Expected results
Results = Array{Float64}(undef, 2, 2, 3)
Results[1,1,1] = -9.1467500
Results[1,2,1] = 5.02825
Results[2,1,1] = 19.20325
Results[2,2,1] = 33.37825
Results[1,1,2] = -8.7930500
Results[1,2,2] = 5.38195
Results[2,1,2] = 19.55695
Results[2,2,2] = 33.73195
Results[1,1,3] = -5.256050
Results[1,2,3] = 8.91895
Results[2,1,3] = 23.09395
Results[2,2,3] = 37.26895
# Creating climgrids
axisdata_tdiu = AxisArray(data_tdiu, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_vp = AxisArray(data_vp, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C_tdiu = ClimateTools.ClimGrid(axisdata_tdiu, dataunits= "K", variable = "tmean")
C_vp = ClimateTools.ClimGrid(axisdata_vp, dataunits = "Pa", variable = "vp")
# Use the function
C_wbgt = wbgt(C_tdiu, C_vp)
# Run the test
@test round.(C_wbgt.data.data, digits=10) == Results

# Test for tdiu in Celsius
data_tdiu = Array{Float64,3}(undef, 2,2,3)
data_tdiu[1,1,:] .= 250.0 - 273.15
data_tdiu[1,2,:] .= 275.0 - 273.15
data_tdiu[2,1,:] .= 300.0 - 273.15
data_tdiu[2,2,:] .= 325.0 - 273.15
# data_tdiu = [data_tdiu][1]u"°C"

axisdata_tdiu = AxisArray(data_tdiu, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C_tdiu = ClimateTools.ClimGrid(axisdata_tdiu, dataunits= "Celsius", variable = "tmean")
# Use the function
C_wbgt = wbgt(C_tdiu, C_vp)
# Run the test
@test round.(C_wbgt.data.data, digits=10) == Results

# Test diurnaltemperature()
d = DateTime(2003,1,1):Day(1):DateTime(2003,1,3)
# Dummy data
data_tmax = Array{Float64,3}(undef, 2,2,3)
data_tmax[1,1,:] .= 0.0
data_tmax[1,2,:] .= 10.0
data_tmax[2,1,:] .= 20.0
data_tmax[2,2,:] .= 30.0
data_tmin = Array{Float64,3}(undef, 2,2,3)
data_tmin[:,:,1] .= -10.0
data_tmin[:,:,2] .= 0.0
data_tmin[:,:,3] .= 10.0
# data_tmin = [data_tmin][1]u"°C"

α = 0.4
# Expected results
Results[1,1,1] = -4.0
Results[1,2,1] = 2.0
Results[2,1,1] = 8.0
Results[2,2,1] = 14.0
Results[1,1,2] = 0.0
Results[1,2,2] = 6.0
Results[2,1,2] = 12.0
Results[2,2,2] = 18.0
Results[1,1,3] = 4.0
Results[1,2,3] = 10.0
Results[2,1,3] = 16.0
Results[2,2,3] = 22.0
# Results = [Results][1]u"°C"

# Creating climgrids
axisdata_tmax = AxisArray(data_tmax, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_tmin = AxisArray(data_tmin, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C_tmax = ClimateTools.ClimGrid(axisdata_tmax, dataunits= "Celsius", variable = "tasmax")
C_tmin = ClimateTools.ClimGrid(axisdata_tmin, dataunits = "Celsius", variable = "tasmin")
# Using the function
C_tdiu = diurnaltemperature(C_tmin, C_tmax, α)
# Run the test
@test C_tdiu.data.data == Results

# Test meantemperature()
d = DateTime(2003,1,1):Day(1):DateTime(2003,1,3)
# Dummy data
data_tmax = Array{Float64,3}(undef, 2,2,3)
data_tmax[1,1,:] .= 0.0
data_tmax[1,2,:] .= 10.0
data_tmax[2,1,:] .= 20.0
data_tmax[2,2,:] .= 30.0
# data_tmax = [data_tmax][1]u"°C"

data_tmin = Array{Float64,3}(undef, 2,2,3)
data_tmin[:,:,1] .= -10.0
data_tmin[:,:,2] .= 0.0
data_tmin[:,:,3] .= 10.0
# data_tmin = [data_tmin][1]u"°C"

# Expected results
Results = Array{Float64,3}(undef, 2,2,3)
Results[1,1,1] = -5.0
Results[1,2,1] = 0.0
Results[2,1,1] = 5.0
Results[2,2,1] = 10.0
Results[1,1,2] = 0.0
Results[1,2,2] = 5.0
Results[2,1,2] = 10.0
Results[2,2,2] = 15.0
Results[1,1,3] = 5.0
Results[1,2,3] = 10.0
Results[2,1,3] = 15.0
Results[2,2,3] = 20.0
# Results = [Results][1]u"°C"

# Creating climgrids
axisdata_tmax = AxisArray(data_tmax, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_tmin = AxisArray(data_tmin, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C_tmax = ClimateTools.ClimGrid(axisdata_tmax, dataunits= "Celsius", variable = "tasmax", frequency="day")
C_tmin = ClimateTools.ClimGrid(axisdata_tmin, dataunits = "Celsius", variable = "tasmin")
# Using the function
C_mean = meantemperature(C_tmin, C_tmax)
# Run the test
@test C_mean.data.data == Results

C_tmax2 = C_tmax
ens = [C_tmax2, C_tmax, C_tmax, C_tmax, C_tmax]
E = ensemble_mean(ens)
@test E[1][1, 1, 1] == C_tmax[1][1, 1, 1]

E = ensemble_std(ens)
@test E[1][1, 1, 1] == std([C_tmax[1][1,1,1], C_tmax[1][1,1,1], C_tmax[1][1,1,1], C_tmax[1][1, 1, 1], C_tmax[1][1,1,1]])

d = DateTime(2003,1,1):Day(1):DateTime(2003,1,3)
# Dummy data
data_tmax = Array{Float64,3}(undef, 2,2,3)
data_tmax[1,1,:] .= 0.0
data_tmax[1,2,:] .= 10.0
data_tmax[2,1,:] .= 20.0
data_tmax[2,2,:] .= 30.0
# data_tmax = [data_tmax][1]u"K"

axisdata_tmax = AxisArray(data_tmax, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C_tmax = ClimateTools.ClimGrid(axisdata_tmax, dataunits= "Kelvin", variable = "tasmax", frequency="day")

ens = [C_tmax2 + 2.0, C_tmax, C_tmax, C_tmax, C_tmax]
E = ensemble_max(ens)
@test E[1][1,1,1] == ens[1][1][1,1,1]

E = ensemble_min(ens)
@test E[1][1,1,1] == ens[2][1][1,1,1]

# Drought code
# core function
temp = 17.0
precip = 0.0
month = 4
dc0 = 15.0
@test drought_dc(precip, temp, dc0, month) == 19.014

dc0 = 19.014
temp = 20.0
precip = 2.4
@test round(drought_dc(precip, temp, dc0, month), digits=3) == 23.568

dc0 = 23.568
temp = 8.5
precip = 0.0
@test round(drought_dc(precip, temp, dc0, month), digits=3) == 26.052

dc0 = 40.6
temp = 7.0
precip = 9.0
@test round(drought_dc(precip, temp, dc0, month), digits=3) == 29.529

# Drought code
# ClimGrid wrapper
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
datevec = collect(DateTime(2000, 4, 13):Day(1):DateTime(2000, 4, 27))
temp = Array{Float32}(undef, 1, 1, 15)
temp[1,1,:] = [17.0, 20.0, 8.5, 6.5, 13.0, 6.0, 5.5, 8.5, 9.5, 7.0, 6.5, 6.0, 13.0, 15.5, 23.0]
precip = Array{Float32}(undef, 1, 1, 15)
precip[1,1,:] = [0.0, 2.4, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 9.0, 1.0, 0.0, 0.0, 0.0, 0.0]

tempax = AxisArray(temp, Axis{:lon}([-88.0]), Axis{:lat}([45.0]), Axis{:time}(datevec))
precax = AxisArray(precip, Axis{:lon}([-88.0]), Axis{:lat}([45.0]), Axis{:time}(datevec))

tas = ClimGrid(tempax, variable="tas", dimension_dict=dimension_dict, dataunits="Celsius")
pr = ClimGrid(precax, variable="pr", dimension_dict=dimension_dict, dataunits="mm")

dc = drought_dc(pr, tas)

res = round.([19.014, 23.567999999999998, 26.052, 28.176, 31.47, 33.504, 35.448,  37.932, 40.596000000000004, 29.52469510340648, 31.648695103406478, 33.68269510340648, NaN, NaN, NaN], digits=3)

@test sum(round.(dc[1][1,1,1:end-3], digits=3) .== res[1:end-3]) == length(res[1:end-3])

@test sum(isnan.(dc[1][1,1,:])) == 3

end
