
# Shapefile test
filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")

polyshp = read(filename,Shapefile.Handle)
x, y = shapefile_coords(polyshp.shapes[1])
P = [x y]
P = P'
lat = NetCDF.ncread(filenc, "lat")
lon = NetCDF.ncread(filenc, "lon")
msk = inpolygrid(C.longrid, C.latgrid, P)
C = nc2julia(filenc, "tas", poly = P)
status, figh = mapclimgrid(C); @test status == true; PyPlot.close()
C = nc2julia(filenc, "tas")
status, figh = mapclimgrid(C, mask = msk); @test status == true; PyPlot.close()
@test length(x) == 6
@test length(y) == 6
@test isnan(x[1])
@test isnan(y[1])

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
C = nc2julia(filename, "tas")
status, figh = mapclimgrid(C);@test status == true;PyPlot.close()
status, figh = mapclimgrid(C, region = "World");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Canada");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "Quebec");@test status == true; PyPlot.close()
status, figh = mapclimgrid(C, region = "NorthAmerica");@test status == true; PyPlot.close()
status, figh = mapclimgrid(annualmax(C), region = "Europe");@test status == true; PyPlot.close()

# precip
C = nc2julia(filename, "pr")
status, figh = mapclimgrid(C);@test status == true; PyPlot.close()
status, figh = mapclimgrid(prcp1(C), region = "World");@test status == true; PyPlot.close()
status, figh = mapclimgrid(prcp1(C), region = "Canada");@test status == true; PyPlot.close()
status, figh = mapclimgrid(prcp1(C), region = "Europe");@test status == true; PyPlot.close()
status, figh = mapclimgrid(prcp1(C), region = "NorthAmerica");@test status == true; PyPlot.close()

# ua wind
filename = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
polyshp = read(filename,Shapefile.Handle)
x, y = shapefile_coords(polyshp.shapes[1])
P = [x y]
P = P'
C = nc2julia(filenc, "ua")
status, figh = mapclimgrid(C, level = 3);@test status == true; PyPlot.close() # given level
status, figh = mapclimgrid(C);@test status == true; PyPlot.close() # feeding a 4D field
status, figh = mapclimgrid(C, poly = P);@test status == true; PyPlot.close() # feeding a 4D field with a polygon
status, figh = mapclimgrid(C, mask = msk);@test status == true; PyPlot.close() # feeding a 4D field with a mask
C = nc2julia(filenc, "ua", poly = P)
status, figh = mapclimgrid(C);@test status == true; PyPlot.close() # feeding a 4D field
