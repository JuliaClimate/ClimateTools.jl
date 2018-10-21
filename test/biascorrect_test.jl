# replstr(x) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80)), MIME("text/plain"), x), x)
# showstr(x) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80)), x), x)

d = Date(1961,1,1):Day(1):Date(1990,12,31)
Random.seed!(42)
data = randn(2, 2, 10957)
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs * 1.05
fut = obs * 1.05

D = qqmap(obs, ref, fut, method = "Additive", detrend=false)
@test round(D[1][1, 1, 1], digits=5) == -0.55603#68761463861
D = qqmap(obs, ref, fut, method = "Additive", detrend=true)
@test round(D[1][1, 1, 1], digits=2) == -0.56#5603#74620422795

Random.seed!(42)
data = randn(2, 2, 10957)
axisdata = AxisArray(data .- minimum(data), Axis{:lon}(1:2), Axis{:lat}(1:2),Axis{:time}(d))
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs * 1.05
fut = obs * 1.05

D = qqmap(obs, ref, fut, method = "Multiplicative", detrend=false)
@test round(D[1][1, 1, 1], digits=5)== 3.78144#13533272675

# dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
# d = Date(1961,1,1):Day(1):Date(1990,12,31)
# Random.seed!(42)
# data = randn(6, 6, 10957)
# axisdata = AxisArray(data, Axis{:lon}(1:6), Axis{:lat}(1:6), Axis{:time}(d))
# obs = ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
# ref = obs - 3

# ITP = qqmaptf(obs, ref, partition=0.5, detrend = false)
# @test round(ITP.itp[Random.rand(1:365)][Random.randn(1)][1], digits=10) == 3.0

# interface
# @test replstr(ITP) == "TransferFunction type with fields *itp*, *method* and *detrend*\nInterpolation array: (365,) transfer functions\nMethod: Additive\nDetrended: false"

# @test showstr(ITP) == "ClimateTools.TransferFunction(Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Interpolations.Gridded{Interpolations.Linear},Tuple{Array{Float64,1}},0},Interpolations.Gridded{Interpolations.Linear},Interpolations.OnGrid,Interpolations.Flat}[[3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]  …  [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0], [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0  …  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]], \"Additive\", false)"

# fut = obs - 3
# D = qqmap(fut, ITP)
# @test D[1][1, 1, 1] == -0.5560268761463862

#
# dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
# d = Date(1961,1,1):Day(1):Date(1990,12,31)
# Random.seed!(42)
# data = randn(6, 6, 10957)
# axisdata = AxisArray(data, Axis{:lon}(1:6), Axis{:lat}(1:6), Axis{:time}(d))
# obs = ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
# ref = obs - 3

# ITP = qqmaptf(obs, ref, partition=0.5, detrend = true)
# @test round(ITP.itp[rand(1:365)][randn(1)][1], digits=1) == 3.0

# fut = obs - 3
# D = qqmap(fut, ITP)
# @test round(D[1][1, 1, 1], digits=2) == -0.56#5599#11004056662

# ITP = qqmaptf(obs, ref, partition=0.5, detrend = true, method="multiplicative")
# @test round(ITP.itp[rand(1:365)][randn(1)][1], digits=1) == 0.0
# fut = obs - 3
# D = qqmap(fut, ITP)
# @test round(D[1][1, 1, 1], digits=2) == -0.56#603#75611239462


# Create a ClimGrid with a clear trend
x = 1:10957
x = x
d = Date(1961,1,1):Day(1):Date(1990,12,31)
data = zeros(2, 2, 10957)
data[1,1,:] .= 0.0 .+ x
data[1,2,:] .= 0.0 .+ 0.2.*x .+ x.^2
data[2,1,:] .= 0.0 .+ 0.04.*x - 0.01 .* x.^3
data[2,2,:] .= 0.0 .+ 0.14.*x.^2 .- 0.032 .* x .^4
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
C = ClimGrid(axisdata)
poly = polyfit(C)
val = polyval(C, poly)
D = C - val
@test D[1] == (C - val)[1]

