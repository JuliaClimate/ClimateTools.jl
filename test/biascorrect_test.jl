d = Date(1961,1,1):Date(1990,12,31)
srand(42)
data = randn(2, 2, 10957)
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs * 1.05
fut = obs * 1.05

D = qqmap(obs, ref, fut, method = "Additive", detrend=false)
@test D[1][1, 1, 1] == -0.5560268761463861
D = qqmap(obs, ref, fut, method = "Additive", detrend=true)
@test D[1][1, 1, 1] == -0.5560274620422795

srand(42)
data = randn(2, 2, 10957)
axisdata = AxisArray(data-minimum(data), Axis{:lon}(1:2), Axis{:lat}(1:2),Axis{:time}(d))
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs * 1.05
fut = obs * 1.05

D = qqmap(obs, ref, fut, method = "Multiplicative", detrend=false)
@test D[1][1, 1, 1] == 3.7814413533272675

dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
d = Date(1961,1,1):Date(1990,12,31)
srand(42)
data = randn(6, 6, 10957)
axisdata = AxisArray(data, Axis{:lon}(1:6), Axis{:lat}(1:6), Axis{:time}(d))
obs = ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs - 3

ITP = qqmaptf(obs, ref, partition=0.5, detrend = false)
@test round(ITP.itp[rand(1:365)][randn(1)][1],10) == 3.0

fut = obs - 3
D = qqmap(fut, ITP)
@test D[1][1, 1, 1] == -0.5560268761463862

#
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
d = Date(1961,1,1):Date(1990,12,31)
srand(42)
data = randn(6, 6, 10957)
axisdata = AxisArray(data, Axis{:lon}(1:6), Axis{:lat}(1:6), Axis{:time}(d))
obs = ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs - 3

ITP = qqmaptf(obs, ref, partition=0.5, detrend = true)
@test round(ITP.itp[rand(1:365)][randn(1)][1],1) == 3.0

fut = obs - 3
D = qqmap(fut, ITP)
@test D[1][1, 1, 1] == -0.5559911004056662


# Create a ClimGrid with a clear trend
x = 1:10957
x = x
d = Date(1961,1,1):Date(1990,12,31)
data = zeros(2, 2, 10957)
data[1,1,:] = 0.0 + x
data[1,2,:] = 0.0 + 0.2.*x + x.^2
data[2,1,:] = 0.0 + 0.04.*x - 0.01 .* x.^3
data[2,2,:] = 0.0 + 0.14.*x.^2 - 0.032 .* x .^4
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimGrid(axisdata)
poly = polyfit(C)
val = polyval(C, poly)
D = C - val
