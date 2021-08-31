# replstr(x) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80)), MIME("text/plain"), x), x)
# showstr(x) = sprint((io,x) -> show(IOContext(io, :limit => true, :displaysize => (24, 80)), x), x)

@testset "Bias correction" begin

d = DateTime(1961,1,1):Day(1):DateTime(1990,12,31)

filedata = joinpath(dirname(@__FILE__), "data", "data_test.h5")
data = h5read(filedata, "random/data")
#Random.seed!(42)
#data = randn(2, 2, 10957)
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs * 1.05
fut = obs * 1.05

D = qqmap(obs, ref, fut, method = "Additive", detrend=false)
@test round(D[1][1, 1, 1], digits=5) == -0.55603#68761463861
D = qqmap(obs, ref, fut, method = "Additive", detrend=true)
@test round(D[1][1, 1, 1], digits=2) == -0.56#5603#74620422795

# ===================================
# Ensure we have similar statistics
# ===================================
d = DateTime(1961,1,1):Day(1):DateTime(1990,12,31)
# filedata = joinpath(dirname(@__FILE__), "data", "data_test.h5")
# data = h5read(filedata, "random/data")
#Random.seed!(42)
#data = randn(2, 2, 10957)
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs + 2.0
fut = obs + 2.0

D = qqmap(obs, ref, fut, method="Additive", detrend=false)
@test mean(obs[1][1,1,:]) - mean(fut[1][1,1,:]) ≈ -1.99999999999
@test mean(obs[1][1,1,:]) - mean(D[1][1,1,:]) ≈ -0.00029395108

ref = obs * 2.0
fut = obs * 2.0
D = qqmap(obs, ref, fut, method="Multiplicative", detrend=false)
@test std(obs[1][1,1,:]) - std(fut[1][1,1,:]) ≈ -1.012096874
@test round(std(obs[1][1,1,:]) - std(D[1][1,1,:]), digits=3) ≈ -0.0#5.233448404e-5

# ============================
# filedata = joinpath(dirname(@__FILE__), "data", "data_test.h5")
# data = h5read(filedata, "random/data")

axisdata = AxisArray(data .- minimum(data), Axis{:lon}(1:2), Axis{:lat}(1:2),Axis{:time}(d))
obs = ClimateTools.ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)
ref = obs * 1.05
fut = obs * 1.05

D = qqmap(obs, ref, fut, method = "Multiplicative", detrend=false)
@test round(D[1][1, 1, 1], digits=2) == 3.78

# Bias extremes
d = DateTime(1961,1,1):Day(1):DateTime(1990,12,31)
lat = 1.0:2.0
lon = 1.0:2.0

# filedata = joinpath(dirname(@__FILE__), "data", "data_test.h5")
# data = h5read(filedata, "random/data")
#Random.seed!(42)
#data = randn(2, 2, 10957)

longrid, latgrid = ndgrid(lon, lat)
axisdata = AxisArray(data, Axis{:lon}(lon), Axis{:lat}(lat), Axis{:time}(d))
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
varattribs = Dict(["standard_name" => "random noise"])
obs = ClimateTools.ClimGrid(axisdata, variable = "psl", typeofvar="psl", longrid=longrid, latgrid=latgrid, dimension_dict=dimension_dict, varattribs=varattribs)
ref = obs * 1.1
fut = obs * 1.05

mu = [36.6208, 37.423]
sigma = [11.555, 10.381]
xi = [0.08393, 0.08393]
gevparams = DataFrame([lat, lon, mu, sigma, xi], [:lat, :lon, :mu, :sigma, :xi])
Dext = biascorrect_extremes(obs, ref, fut, detrend=true, gevparams=gevparams)
#@test round(maximum(Dext), digits=5) == 120.01175
@test round(maximum(Dext), digits=5) == 103.03353

Dext = biascorrect_extremes(obs, ref, fut, detrend=false, gevparams=gevparams)
#@test round(maximum(Dext), digits=5) == 120.00829
@test round(maximum(Dext), digits=5) == 103.03378
Dext = biascorrect_extremes(obs, ref, fut, detrend=false)
#@test round(maximum(Dext), digits=5) == 16.7384
@test round(maximum(Dext), digits=5) == 4.10549

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

# Test for threshold inclusion or not in Extremes (changing behaviors between version)
filegamma = joinpath(dirname(@__FILE__), "data", "gamma.h5")
obs = h5read(filegamma, "obsgamma")

thres = quantile(obs, 0.95)
clusters = Extremes.getcluster(obs, thres, 0.1)

GPD = Extremes.gpdfit(clusters[!,:Max], threshold = thres)
@test round(quantile(GPD, 0.0), digits=5) .== round(thres, digits=5)
@test round(quantile(GPD, 0.999), digits=5) .== round(123.02851, digits=5)


end
