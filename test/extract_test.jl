# Period subset
d = DateTime(1961,1,1):Day(1):DateTime(1990,12,31)
Random.seed!(42)
data = randn(2, 2, 10957)
axisdata = AxisArray(data, Axis{:lon}(1:2), Axis{:lat}(1:2), Axis{:time}(d))
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
C = ClimateTools.ClimGrid(axisdata, variable = "tasmax", dimension_dict=dimension_dict)

# When startmonth < endmonth
D = resample(C, 3, 6)
@test D[1].data[:,:,1] == C[1].data[:,:,60]
@test length(D[1][Axis{:time}]) == 3660 && Dates.month(D[1][Axis{:time}][end]) == 6 && Dates.month(D[1][Axis{:time}][1]) == 3

# When startmonth > endmonth
D = resample(C, 10, 2)
@test D[1].data[:,:,1] == C[1].data[:,:,1] && D[1].data[:,:,60] == C[1].data[:,:,274]
@test length(D[1][Axis{:time}]) == 4537 && Dates.month.(D[1][Axis{:time}][end]) == 12 && Dates.month.(D[1][Axis{:time}][1]) == 1
