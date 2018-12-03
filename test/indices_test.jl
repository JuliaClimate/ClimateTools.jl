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
# @test prcp1(data, d) == Results

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
# @test annualsum(data, d) == Results

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
# @test frostdays(data, d) == Results
# @test icingdays(data, d) == Results

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
# @test summerdays(data, d) == Results

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
# @test tropicalnights(data, d) == Results

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
# @test customthresover(data, d, 20) == Results

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
# @test customthresunder(data, d, 200) == Results

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
# @test annualmax(data, d) == Results

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
# @test annualmin(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tasmin")
ind = annualmin(C)
@test ind.data.data == Results

# Days above 10

data = Array{Float64, 3}(undef, 2, 2, 1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:] = collect(1.0:1826.0); data[2,1,:]=collect(1.0:1826.0); data[2,2,:] = collect(1.0:1826.0);
Results = Array{Int64, 3}(undef, 2, 2, 5);
Results[1, 1, :]=[356, 366, 365, 365, 365]'';
Results[1, 2, :] = [356, 366, 365, 365, 365]'';
Results[2, 1, :] = [356, 366, 365, 365, 365]'';
Results[2, 2, :] = [356, 366, 365, 365, 365]'';
# @test summerdays(data, d) == Results

# ClimGrid based tests
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimateTools.ClimGrid(axisdata, variable = "tas")
ind = daysabove10(C)
@test ind.data.data == Results

# Period mean
data = Array{Float64, 3}(undef, 2, 2, 1826)
data[1,1,:] = collect(1.0:1826.0); data[1,2,:] = collect(-10.0:1815.0); data[2,1,:] = collect(0.0:1825.0); data[2,2,:] = collect(-1725.0:100);
Results = Array{Float64, 2}(undef, 2, 2);
Results[1, 1]= 746.0;
Results[1, 2] = 735.0;
Results[2, 1] = 745.0;
Results[2, 2] = -980.0;

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
ind = periodmean(C)
@test ind.data.data == Results

# Daymean
raw = collect(0:0.25:1459.75) # 6hr data
d = timedecode(raw, "days since 2000-01-01", "noleap")

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
data_orog = Array{Float64,2}(undef, 2,2)
data_orog[:,1] .= 0.0
data_orog[:,2] .= 10.0
data_tas = Array{Float64,3}(undef, 2,2,3)
data_tas[1,:,:] .= 100.0
data_tas[2,:,:] .= 300.0
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
data_vp = Array{Float64,3}(undef, 2,2,3)
data_vp[:,:,1] .= 10.0
data_vp[:,:,2] .= 100.0
data_vp[:,:,3] .= 1000.0
# Expected results
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
data_tmin = Array{Float64,3}(undef, 2,2,3)
data_tmin[:,:,1] .= -10.0
data_tmin[:,:,2] .= 0.0
data_tmin[:,:,3] .= 10.0
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
ens = [C_tmax2, C_tmax]
E = ensemble_mean(ens)
@test E[1][1, 1, 1] == C_tmax[1][1, 1, 1]
ens = [C_tmax2, C_tmax*2]
@test E[1][1, 1, 1] == (C_tmax[1][1, 1, 1]*2 + C_tmax2[1][1, 1, 1])/2


# d = DateTime(2003,1,1):Hour(1):DateTime(2003,1,3)-Hour(1)
# data = Array{Float64}(undef, 2, 2, 48)
# vecdata = collect(range(1, length=48))

# data[1, 1, :] .= vecdata
# data[1, 2, :] .= vecdata
# data[2, 1, :] .= vecdata
# data[2, 2, :] .= vecdata

# axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
# C = ClimGrid(axisdata, dataunits = "Celsius", variable = "tasmin", frequency="1h")

end
