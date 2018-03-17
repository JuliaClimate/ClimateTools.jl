d = Date(1961,1,1):Date(1990,12,31)
srand(42)
data = randn(10957, 2, 2)
axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax")
ref = obs * 1.05
fut = obs * 1.05

D = qqmap(obs, ref, fut, method = "Additive")
@test D[1][1, 1, 1] == -0.5560268761463861

srand(42)
data = randn(10957, 2, 2)
axisdata = AxisArray(data-minimum(data), Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax")
ref = obs * 1.05
fut = obs * 1.05

D = qqmap(obs, ref, fut, method = "Multiplicative")
@test D[1][1, 1, 1] == 3.781441353327268
