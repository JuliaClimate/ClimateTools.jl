d = Date(1961,1,1):Date(1990,12,31)
srand(42)
data = randn(2, 2, 10957)
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs * 1.05
fut = obs * 1.05

D = qqmap(obs, ref, fut, method = "Additive")
@test D[1][1, 1, 1] == -0.5560268761463861

srand(42)
data = randn(2, 2, 10957)
axisdata = AxisArray(data-minimum(data), Axis{:lon}(1:2), Axis{:lat}(1:2),Axis{:time}(d))
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs * 1.05
fut = obs * 1.05

D = qqmap(obs, ref, fut, method = "Multiplicative")
@test D[1][1, 1, 1] == 3.7814413533272675
