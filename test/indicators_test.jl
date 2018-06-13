# Test vaporpressure(specific_humidity::ClimGrid, surface_pressure::ClimGrid)
d = Date(2003,1,1):Date(2003,1,3)
# Dummy data
data_huss = Array{Float64,3}(2,2,3)
data_huss[1,1,:] = 0
data_huss[1,2,:] = 0.1
data_huss[2,1,:] = 0.5
data_huss[2,2,:] = 1
data_ps = Array{Float64,3}(2,2,3)
data_ps[:,:,1] = 10
data_ps[:,:,2] = 100
data_ps[:,:,3] = 1000
# Expected values
Results = Array{Float64,3}(2,2,3)
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
data_huss = Array{Float64,3}(2,2,3)
data_huss[1,1,:] = 0
data_huss[1,2,:] = 0.1
data_huss[2,1,:] = 0.5
data_huss[2,2,:] = 1
data_psl = Array{Float64,3}(2,2,3)
data_psl[:,:,1] = 10
data_psl[:,:,2] = 100
data_psl[:,:,3] = 1000
data_orog = Array{Float64,2}(2,2)
data_orog[:,1] = 0
data_orog[:,2] = 10
data_tas = Array{Float64,3}(2,2,3)
data_tas[1,:,:] = 100
data_tas[2,:,:] = 300
# Expected results
Results = Array{Float64,3}(2,2,3)
Results[1,1,:] = round(0.0, 10)
Results[1,2,1] = round(1.380315267090330, 10)
Results[1,2,2] = round(13.80315267090330, 10)
Results[1,2,3] = round(138.0315267090330, 10)
Results[2,1,1] = round(4.456327985739750, 10)
Results[2,1,2] = round(44.56327985739750, 10)
Results[2,1,3] = round(445.6327985739750, 10)
Results[2,2,1] = round(6.158207427095860, 10)
Results[2,2,2] = round(61.58207427095860, 10)
Results[2,2,3] = round(615.8207427095860, 10)
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
@test round.(vp.data.data, 10) == Results

# Test for tas in Celsius
data_tas = Array{Float64,3}(2,2,3)
data_tas[1,:,:] = 100-273.15
data_tas[2,:,:] = 300-273.15
axisdata_tas = AxisArray(data_tas, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C_tas = ClimateTools.ClimGrid(axisdata_tas, dataunits = "Celsius", variable = "tas")
# Using the function on dummy data
vp = vaporpressure(C_huss, C_psl, C_orog, C_tas)
# Run the test
@test round.(vp.data.data, 10) == Results

# Test wbgt(diurnal_temperature::ClimGrid, vapor_pressure::ClimGrid)
d = Date(2003,1,1):Date(2003,1,3)
# Dummy data
data_tdiu = Array{Float64,3}(2,2,3)
data_tdiu[1,1,:] = 250
data_tdiu[1,2,:] = 275
data_tdiu[2,1,:] = 300
data_tdiu[2,2,:] = 325
data_vp = Array{Float64,3}(2,2,3)
data_vp[:,:,1] = 10
data_vp[:,:,2] = 100
data_vp[:,:,3] = 1000
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
@test round.(C_wbgt.data.data,10) == Results

# Test for tdiu in Celsius
data_tdiu = Array{Float64,3}(2,2,3)
data_tdiu[1,1,:] = 250-273.15
data_tdiu[1,2,:] = 275-273.15
data_tdiu[2,1,:] = 300-273.15
data_tdiu[2,2,:] = 325-273.15
axisdata_tdiu = AxisArray(data_tdiu, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C_tdiu = ClimateTools.ClimGrid(axisdata_tdiu, dataunits= "Celsius", variable = "tdiu")
# Use the function
C_wbgt = wbgt(C_tdiu, C_vp)
# Run the test
@test round.(C_wbgt.data.data,10) == Results

# Test diurnaltemperature()
d = Date(2003,1,1):Date(2003,1,3)
# Dummy data
data_tmax = Array{Float64,3}(2,2,3)
data_tmax[1,1,:] = 0.0
data_tmax[1,2,:] = 10.0
data_tmax[2,1,:] = 20.0
data_tmax[2,2,:] = 30.0
data_tmin = Array{Float64,3}(2,2,3)
data_tmin[:,:,1] = -10.0
data_tmin[:,:,2] = 0.0
data_tmin[:,:,3] = 10.0
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
d = Date(2003,1,1):Date(2003,1,3)
# Dummy data
data_tmax = Array{Float64,3}(2,2,3)
data_tmax[1,1,:] = 0.0
data_tmax[1,2,:] = 10.0
data_tmax[2,1,:] = 20.0
data_tmax[2,2,:] = 30.0
data_tmin = Array{Float64,3}(2,2,3)
data_tmin[:,:,1] = -10.0
data_tmin[:,:,2] = 0.0
data_tmin[:,:,3] = 10.0
# Expected results
Results = Array{Float64,3}(2,2,3)
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
C_tmax = ClimateTools.ClimGrid(axisdata_tmax, dataunits= "Celsius", variable = "tasmax")
C_tmin = ClimateTools.ClimGrid(axisdata_tmin, dataunits = "Celsius", variable = "tasmin")
# Using the function
C_mean = meantemperature(C_tmin, C_tmax)
# Run the test
@test C_mean.data.data == Results
