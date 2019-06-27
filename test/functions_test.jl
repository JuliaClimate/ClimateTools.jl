using Dates, AxisArrays

@testset "Functions" begin
# findmax
d = collect(DateTime(2003,1,1):Day(1):DateTime(2005,12,31))

data = Array{Float64,3}(undef, 2, 2,1096)
data[1,1,:] = collect(1.0:1096.0); data[1,2,:] = collect(1.0:1096.0); data[2,1,:]=collect(1.0:1096.0); data[2,2,:] = collect(1.0:1096.0);
data[1,1,1] = NaN

# ClimGrid based tests
x = 1:2
y = 1:2
axisdata = AxisArray(data, Axis{:lon}(x), Axis{:lat}(y), Axis{:time}(d))
X, Y = ndgrid(x, y)
C = ClimateTools.ClimGrid(axisdata, variable = "pr", longrid=X, latgrid=Y)
@test isnan(findmin(C)[1])
@test findmin(C)[2] == CartesianIndex(1,1,1)
@test isnan(findmax(C)[1])
@test findmax(C)[2] == CartesianIndex(1,1,1)
@test findmin(C, skipnan=true) == (1.0, CartesianIndex(2,1,1))
@test findmax(C, skipnan=true) == (1096.0, CartesianIndex(1,1,1096))

# Monthly means
D = monthmean(C)
@test size(D[1]) == (2, 2, 36)
@test isnan(D[1][1, 1, 1])
@test D[1][2,1,1] == 16.0
@test D[1][2,2,1] == 16.0
@test D[1][2,2,end] == 1081.0
timevec = get_timevec(D)
@test timevec[1] == DateTime(2003,01,01,00,00,00)
@test timevec[end] == DateTime(2005,12,01,00,00,00)

D = monthsum(C)
@test size(D[1]) == (2, 2, 36)
@test isnan(D[1][1, 1, 1])
@test D[1][2,1,1] == 496.0
@test D[1][2,2,1] == 496.0
@test D[1][2,2,end] == 33511.0
timevec = get_timevec(D)
@test timevec[1] == DateTime(2003,01,01,00,00,00)
@test timevec[end] == DateTime(2005,12,01,00,00,00)

# get dims
x, y, timevec = ClimateTools.getdims(C)
@test length(timevec) == length(d)
for itime = 1:length(timevec)
    @test timevec[itime] == d[itime]
end
@test x == 1:2 && y == 1:2

# get grids
longrid, latgrid = ClimateTools.getgrids(C)
@test longrid == X
@test latgrid == Y

# Extension file
@test ClimateTools.extension("test") == ""
@test ClimateTools.extension("test.nc") == ".nc"
@test ClimateTools.extension("test.") == ""

dt = DateTimeNoLeap(2000, 01, 01, 23)
@test yearmonthdayhour(dt) == (2000, 1, 1, 23)

# WeatherStation based tests
d_ws = collect(DateTime(2003,1,1):Year(1):DateTime(2005,12,31))
data_ws = Array{Float64,1}(undef, 3)
axisdata_ws = AxisArray(data_ws, Axis{Symbol(time)}(d2))
W = WeatherStation(axisdata_ws, Real(-65), Real(45), Real(30), "TEST000A")
timevec_ws = ClimateTools.get_timevec(W)
@test timevec_ws[1] == DateTime(2003,01,01,00,00,00)
@test timevec_ws[end] == DateTime(2005,01,01,00,00,00)

end
