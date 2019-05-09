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

maps = ["rotated_pole", "rotated_mercator", "crs", "lambert_conformal_conic","rotated_latitude_longitude", "None"]

for imap in maps
  K = ["time", imap]
  if imap != "None"
    @test ClimateTools.get_mapping(K) == imap
  else
    @test ClimateTools.get_mapping(K) == "Regular_longitude_latitude"
  end
end

end

@testset "Exportation" begin

    filenc = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
    C = load(filenc, "tas")

    status = write(C, "test.nc")
    @test status == true
    run(`rm test.nc`)

    status = write(C, "test")
    @test status == true
    run(`rm test.nc`)

    # Custom ClimGrid
    data = Float32.(rand(2,2, 365))
    timev = DateTimeNoLeap(2001, 01, 01):Day(1):DateTimeNoLeap(2001, 12, 31)
    rlon = collect(range(20.0, stop=21.0))
    rlat = collect(range(40.0, stop=41.0))
    longrid, latgrid = ndgrid(rlon, rlat)
    variable = "pr"
    A = AxisArray(data, Axis{:rlon}(rlon), Axis{:rlat}(rlat), Axis{:time}(timev))
    dimension_dict = Dict(["lat" => "rlat", "lon" => "rlon"])

    timeattrib = Dict(["units" => "days since 2001-1-1",
    "calendar" => "365_day",
    "axis" => "T",
    "standard_name" => "time"])

    varattrib = Dict(["standard_name" => "precipitation_flux",
    "units" => "kg m-2 s-1"])

    globalattribs = Dict(["title"                 => "CanESM2 model output prepared for CMIP5 historical",
  "branch_time"           => 135415.0,
  "physics_version"       => 1,
  "tracking_id"           => "97542c1e-3841-4a49-bc21-76d6edb9fa1d",
  "branch_time_YMDH"      => "2221:01:01:00",
  "modeling_realm"        => "atmos",
  "creation_date"         => "2011-04-09T02:32:12Z",
  "history"               => "",
  "table_id"              => "Table day (28 March 2011) f9d6cfec5981bb8be1801b35a81002f0",
  "CCCma_runid"           => "IDE",
  "project_id"            => "CMIP5",
  "experiment"            => "historical",
  "forcing"               => "GHG,Oz,SA,BC,OC,LU,Sl,Vl (GHG includes CO2,CH4,N2O,CFC11,effective CFC12)",
  "contact"               => "cccma_info@ec.gc.ca",
  "Conventions"           => "CF-1.4",
  "institute_id"          => "CCCma",
  "parent_experiment_rip" => "r1i1p1",
  "experiment_id"         => "historical",
  "model_id"              => "CanESM2",
  "frequency"             => "day"])

  grid_mapping = Dict(["grid_mapping" => "Rotated_pole"])

  C = ClimGrid(A, variable=variable, longrid=longrid, latgrid=latgrid, dimension_dict=dimension_dict, timeattrib=timeattrib, varattribs=varattrib, globalattribs=globalattribs, grid_mapping=grid_mapping)

  status = write(C, "test.nc")
  @test status == true
  run(`rm test.nc`)

  # grid_mapping = Dict(["grid_mapping_name" => "Rotated_pole"])
  #
  # C = ClimGrid(A, variable=variable, longrid=longrid, latgrid=latgrid, dimension_dict=dimension_dict, timeattrib=timeattrib, varattribs=varattrib, globalattribs=globalattribs, grid_mapping=grid_mapping)
  #
  # status = write(C, "test.nc")
  # @test status == true
  # run(`rm test.nc`)

end
