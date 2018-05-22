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
