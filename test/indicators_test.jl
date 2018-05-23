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
C_ps = ClimateTools.ClimGrid(axisdata_ps, variable = "ps")
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
# Using the function on dummy data
axisdata_huss = AxisArray(data_huss, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_psl = AxisArray(data_psl, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_tas = AxisArray(data_tas, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_psl = AxisArray(data_psl, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
axisdata_orog = AxisArray(data_orog, Axis{:lon}(1:2), Axis{:lat}(1:2))
C_huss = ClimateTools.ClimGrid(axisdata_huss, variable = "huss")
C_psl = ClimateTools.ClimGrid(axisdata_psl, variable = "psl")
C_tas = ClimateTools.ClimGrid(axisdata_tas, variable = "tas")
C_orog = ClimateTools.ClimGrid(axisdata_orog, variable = "orog")
vp = vaporpressure(C_huss, C_psl, C_orog, C_tas)
# Run the test
@test round.(vp.data.data, 10) == Results
