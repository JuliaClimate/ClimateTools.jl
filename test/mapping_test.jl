@testset "Mapping" begin

# Shapefile test
filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")

polyshp = read(filename,Shapefile.Handle)
P = shapefile_coords_poly(polyshp.shapes[1])
# P = [x y]
# P = P'
lat = NetCDF.ncread(filenc, "lat")
lon = NetCDF.ncread(filenc, "lon")
C = load(filenc, "tas", poly = P)
msk = inpolygrid(C.longrid, C.latgrid, P)
status, figh = mapclimgrid(C); @test status == true; PyPlot.close()
status, figh = mapclimgrid(C, mask = msk); @test status == true; PyPlot.close()
C = load(filenc, "tas")
status, figh = mapclimgrid(C); @test status == true; PyPlot.close()
@test size(P) == (2, 6)
@test isnan(P[1,1])
@test isnan(P[2,1])

# filename = joinpath(dirname(@__FILE__), "data", "zoneAgricoleQc15km.shp")
# polyshp = read(filename,Shapefile.Handle)
# x, y = shapefile_coords(polyshp.shapes[2])
# P = [x y]
# P = P'
# cities_lat = collect(45.5016889:0.1:46.8138783)
# cities_lon = collect(-73.56725599999999:0.1:-71.2079809) + 360
# @test sum(.!isnan.(inpolyvec(cities_lon, cities_lat, P))) == 336
# @test sum(.!isnan.(inpolyvec(cities_lon - 360, cities_lat, P))) == 336

# Mapping test
filename = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
C = load(filename, "tas")
status, figh = mapclimgrid(C);@test status == true;PyPlot.close()
status, figh = mapclimgrid(C, region = "World");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "WorldAz");@test status == true; PyPlot.close()
# status, figh = mapclimgrid(C, region = "WorldEck4");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Canada");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Quebec");@test status == true; PyPlot.close()
# status, figh = mapclimgrid(C, region = "QuebecNSP");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Americas");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "NorthAmerica");@test status == true; PyPlot.close()
# status, figh = mapclimgrid(C, region = "SouthAmerica");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Greenwich");@test status == true; PyPlot.close()
status, figh = mapclimgrid(annualmax(C), region = "Europe");@test status == true; PyPlot.close()
# status, figh = mapclimgrid(annualmin(C), region = "Arctic");@test status == true; PyPlot.close()
# status, figh = mapclimgrid(annualmin(C), region = "Antarctic");@test status == true; PyPlot.close()

# precip
C = load(filename, "pr", data_units="mm") + 2.0
status, figh = mapclimgrid(C);@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, start_date=(2000, 05, 16, 12), end_date=(2000, 05, 16, 12));@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "World");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "WorldAz");@test status == true; PyPlot.close()
# status, figh = mapclimgrid(C, region = "WorldEck4");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Canada");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Quebec");@test status == true; PyPlot.close()
# status, figh = mapclimgrid(C, region = "QuebecNSP");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Americas");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "NorthAmerica");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "SouthAmerica");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Americas");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Arctic");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Antarctic");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Greenwich");@test status == true; PyPlot.close()
status, figh = mapclimgrid(annualmax(C), region = "Europe");@test status == true; PyPlot.close()
# status, figh = mapclimgrid(annualmin(C), region = "Arctic");@test status == true; PyPlot.close()
# status, figh = mapclimgrid(annualmin(C), region = "Antarctic");@test status == true; PyPlot.close()

# ua wind
filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
polyshp = read(filename,Shapefile.Handle)
x, y = shapefile_coords(polyshp.shapes[1])
P = [x y]
P = P'
C = load(filenc, "ua")
status, figh = mapclimgrid(C, level = 3);@test status == true; PyPlot.close() # given level
status, figh = mapclimgrid(C);@test status == true; PyPlot.close() # feeding a 4D field
status, figh = mapclimgrid(C, poly = P);@test status == true; PyPlot.close() # feeding a 4D field with a polygon
status, figh = mapclimgrid(spatialsubset(C, P), mask = msk);@test status == true; PyPlot.close() # feeding a 4D field with a mask
C = load(filenc, "ua", poly = P)
status, figh = mapclimgrid(C);@test status == true; PyPlot.close() # feeding a 4D field

regions = ["World", "WorldAz", "Canada", "Quebec", "Quebec_agricole", "South_Quebec", "QuebecNSP", "Americas", "NorthAmerica", "SouthAmerica", "Greenwich", "Outaouais", "Laurentides", "Estrie", "Arctic", "Antarctic", "Africa", "Europe"]

for iregion in regions
    println(iregion)
    status, figh = mapclimgrid(region=iregion);@test status == true; PyPlot.close()
end

# DUMMY MAPS
lon = collect(-180.0:180.0)
lat = collect(-90.0:90.0)
longrid, latgrid = ndgrid(lon, lat)
timeV = DateTime(2003,1,1):Day(1):DateTime(2003,01,31)
data = randn(361, 181, 31)
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
varattribs = Dict(["standard_name" => "random noise"])
axisdata = AxisArray(data, Axis{:lon}(lon), Axis{:lat}(lat), Axis{:time}(timeV))
C = ClimateTools.ClimGrid(axisdata, variable = "psl", typeofvar="psl", longrid=longrid, latgrid=latgrid, dimension_dict=dimension_dict, varattribs=varattribs)

status, figh = mapclimgrid(C);@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, center_cs=true, caxis=[200, 250]);@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, center_cs=true, surface=:pcolormesh);@test status == true; PyPlot.close()


lon = collect(-180.0:180.0)
lat = collect(-90.0:90.0)
longrid, latgrid = ndgrid(lon, lat)
timeV = Dates.year(DateTime(2003, 01, 01)):Dates.year(1):Dates.year(DateTime(2004, 01, 01))
data = randn(361, 181, 2)
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
varattribs = Dict(["standard_name" => "random noise"])
axisdata = AxisArray(data, Axis{:lon}(lon), Axis{:lat}(lat), Axis{:time}(timeV))
C = ClimateTools.ClimGrid(axisdata, variable = "psl", typeofvar="psl", longrid=longrid, latgrid=latgrid, dimension_dict=dimension_dict, varattribs=varattribs)

status, figh = mapclimgrid(C);@test status == true; PyPlot.close()

lon = collect(-180.0:180.0)
lat = collect(-90.0:90.0)
longrid, latgrid = ndgrid(lon, lat)
timeV = DateTime(2003, 01, 01):Day(1):DateTime(2004, 01, 01)
data = randn(361, 181, 366)
dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
varattribs = Dict(["standard_name" => "random noise"])
axisdata = AxisArray(data, Axis{:lon}(lon), Axis{:lat}(lat), Axis{:time}(timeV))
C = ClimateTools.ClimGrid(axisdata, variable = "psl", typeofvar="psl", longrid=longrid, latgrid=latgrid, dimension_dict=dimension_dict, varattribs=varattribs)

status, figh = plot(C); @test status == true; PyPlot.close()
status, figh = plot(C, poly=P); @test status == true; PyPlot.close()
status, figh = plot(C, start_date=(2003, 01, 10), end_date=(2003, 05, 12)); @test status==true;PyPlot.close()
status, figh = plot(C, xlimit=[2 5]); @test status==true;PyPlot.close()
status, figh = plot(C, label = "dummy", titlestr="dummytest", gridfig=true); @test status == true; PyPlot.close()

status, figh = hist(C); @test status == true; PyPlot.close()
status, figh = hist(C, poly=P, ylimit=[2, 5]); @test status == true; PyPlot.close()
status, figh = hist(C, start_date=(2003, 01, 10), end_date=(2003, 05, 12)); @test status==true;PyPlot.close()
status, figh = hist(C, label = "dummy", titlestr="dummytest", gridfig=false); @test status == true; PyPlot.close()

end
