@testset "Extraction" begin
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
@test unique(Dates.month.(D[1][Axis{:time}][:])) == [1,2,10,11,12]

D = resample(C, "DJF") # hardcoded seasons
@test D[1].data[:,:,60] == C[1].data[:,:,365-30]
@test length(D[1][Axis{:time}]) == 2707 && Dates.month(D[1][Axis{:time}][1]) == 1 && Dates.month(D[1][Axis{:time}][end]) == 12
@test unique(Dates.month.(D[1][Axis{:time}][:])) == [1,2,12]

D = resample(C, "MAM") # hardcoded seasons
@test D[1].data[:,:,1] == C[1].data[:,:,60]
@test length(D[1][Axis{:time}]) == 2760 && Dates.month(D[1][Axis{:time}][1]) == 3 && Dates.month(D[1][Axis{:time}][end]) == 5
@test unique(Dates.month.(D[1][Axis{:time}][:])) == [3,4,5]

D = resample(C, "JJA") # hardcoded seasons
@test D[1].data[:,:,1] == C[1].data[:,:,152]
@test length(D[1][Axis{:time}]) == 2760 && Dates.month(D[1][Axis{:time}][1]) == 6 && Dates.month(D[1][Axis{:time}][end]) == 8
@test unique(Dates.month.(D[1][Axis{:time}][:])) == [6,7,8]

D = resample(C, "SON") # hardcoded seasons
@test D[1].data[:,:,1] == C[1].data[:,:,244]
@test length(D[1][Axis{:time}]) == 2730 && Dates.month(D[1][Axis{:time}][1]) == 9 && Dates.month(D[1][Axis{:time}][end]) == 11
@test unique(Dates.month.(D[1][Axis{:time}][:])) == [9,10,11]

@test_throws ErrorException D = resample(C, "df")

end
