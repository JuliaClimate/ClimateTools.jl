@testset "Bias correction tests" begin
# replstr(x) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80)), MIME("text/plain"), x), x)
# showstr(x) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80)), x), x)

@testset "Bias-correction" begin

d = DateTime(1961,1,1):Day(1):DateTime(1990,12,31)
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

# Create a ClimGrid with a clear trend
x = 1:10957
x = x
d = DateTime(1961,1,1):Day(1):DateTime(1990,12,31)
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
