@testset "spatialsubset" begin
    lon = collect(-4.0:1.0:4.0)
    lat = collect(41.0:1.0:45.0)
    time = collect(DateTime(2001, 1, 1):Day(1):DateTime(2001, 1, 2))

    data = reshape(Float32.(1:length(lon) * length(lat) * length(time)), length(lon), length(lat), length(time))
    cube = YAXArray((Dim{:longitude}(lon), Dim{:latitude}(lat), Dim{:time}(time)), data)

    # Triangle so that the bounding-box subset contains both valid and masked points.
    poly = [NaN -2.0 2.0 0.0 -2.0; NaN 41.0 41.0 44.0 41.0]

    sub = spatialsubset(cube, poly)
    @test size(sub) == (4, 3, 2)

    subarr = Array(sub)
    @test any(isnan, subarr[:, :, 1])
    @test any(isfinite, subarr[:, :, 1])

    sub_t = spatialsubset(cube, poly')
    @test size(sub_t) == size(sub)

    arr_t = Array(sub_t)
    @test all((isnan.(subarr) .& isnan.(arr_t)) .| (subarr .== arr_t))

    legacy_like = Any[cube, "K", "N/A"]
    legacy_sub = spatialsubset(legacy_like, poly)
    @test legacy_sub[2] == "K"
    @test size(legacy_sub[1]) == (4, 3, 2)

    # Real shapefile polygon path (SudQC) on YAXArrays cube.
    shpfile = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
    shp_poly = ClimateTools.extractpoly(shpfile, n=1)

    lon_qc = collect(-90.0:1.0:-50.0)
    lat_qc = collect(40.0:1.0:55.0)
    time_qc = collect(DateTime(2001, 1, 1):Day(1):DateTime(2001, 1, 3))
    data_qc = reshape(Float32.(1:length(lon_qc) * length(lat_qc) * length(time_qc)), length(lon_qc), length(lat_qc), length(time_qc))
    cube_qc = YAXArray((Dim{:longitude}(lon_qc), Dim{:latitude}(lat_qc), Dim{:time}(time_qc)), data_qc)

    sub_qc = spatialsubset(cube_qc, shp_poly)
    @test size(sub_qc, 1) < size(cube_qc, 1)
    @test size(sub_qc, 2) <= size(cube_qc, 2)

    arr_qc = Array(sub_qc)
    @test any(isfinite, arr_qc[:, :, 1])
    @test all(isfinite, arr_qc)

    sub_qc_t = spatialsubset(cube_qc, shp_poly')
    arr_qc_t = Array(sub_qc_t)
    @test size(sub_qc_t) == size(sub_qc)
    @test all((isnan.(arr_qc) .& isnan.(arr_qc_t)) .| (arr_qc .== arr_qc_t))

    # Same region on 0-360 longitude coordinates to verify reference alignment.
    lon_qc_360 = map(x -> x < 0 ? x + 360 : x, lon_qc)
    cube_qc_360 = YAXArray((Dim{:longitude}(lon_qc_360), Dim{:latitude}(lat_qc), Dim{:time}(time_qc)), data_qc)
    sub_qc_360 = spatialsubset(cube_qc_360, shp_poly)

    @test size(sub_qc_360) == size(sub_qc)
    arr_qc_360 = Array(sub_qc_360)
    @test all((isnan.(arr_qc) .& isnan.(arr_qc_360)) .| (arr_qc .== arr_qc_360))

    # Disk-backed cube should stay lazy after spatialsubset.
    ncfile = joinpath(dirname(@__FILE__), "data", "sresa1b_ncar_ccsm3-example.nc")
    cube_disk = Cube(open_dataset(ncfile))
    sub_disk = spatialsubset(cube_disk, shp_poly)
    @test !(sub_disk.data isa Array)
    @test occursin("BroadcastDiskArray", string(typeof(sub_disk.data)))

    north_pole_lon = 0.0
    north_pole_lat = 42.5

    rlon_vals = collect(0.0:1.0:5.0)
    rlat_vals = collect(-2.0:1.0:3.0)
    rlon2d, rlat2d = ClimateTools.ndgrid(rlon_vals, rlat_vals)
    lon2d, lat2d = ClimateTools.rotated_to_geographic(rlon2d, rlat2d, north_pole_lon, north_pole_lat)

    rotated_data = reshape(Float32.(1:length(rlon_vals) * length(rlat_vals)), length(rlon_vals), length(rlat_vals))
    rotated_cube = YAXArray(
        (Dim{:rlon}(rlon_vals), Dim{:rlat}(rlat_vals)),
        rotated_data,
        Dict{String, Any}("grid_mapping" => "rotated_pole"),
    )
    lon_cube = YAXArray((Dim{:rlon}(rlon_vals), Dim{:rlat}(rlat_vals)), Float64.(lon2d))
    lat_cube = YAXArray((Dim{:rlon}(rlon_vals), Dim{:rlat}(rlat_vals)), Float64.(lat2d))
    rp_cube = YAXArray(
        (Dim{:maxStrlen64}(1:1),),
        [' '],
        Dict{String, Any}(
            "grid_north_pole_longitude" => north_pole_lon,
            "grid_north_pole_latitude" => north_pole_lat,
        ),
    )
    rotated_ds = Dataset(pr=rotated_cube, lon=lon_cube, lat=lat_cube, rotated_pole=rp_cube)

    rotated_poly = [NaN 1.0 4.0 4.0 1.0 1.0; NaN -1.0 -1.0 2.0 2.0 -1.0]
    geo_poly = copy(rotated_poly)
    valid_vertices = .!isnan.(rotated_poly[1, :])
    vertex_pairs = [
        begin
            lon_vertex, lat_vertex = ClimateTools.rotated_to_geographic(
                reshape([rotated_poly[1, index]], 1, 1),
                reshape([rotated_poly[2, index]], 1, 1),
                north_pole_lon,
                north_pole_lat,
            )
            (lon_vertex[1, 1], lat_vertex[1, 1])
        end for index in findall(valid_vertices)
    ]
    geo_poly[1, valid_vertices] = first.(vertex_pairs)
    geo_poly[2, valid_vertices] = last.(vertex_pairs)

    rotated_sub = spatialsubset(rotated_ds, geo_poly)

    expected_mask = ClimateTools.inpolygrid(Float64.(lon2d), Float64.(lat2d), geo_poly)
    expected_points = findall(!isnan, expected_mask)
    expected_lon_range = minimum(getindex.(expected_points, 1)):maximum(getindex.(expected_points, 1))
    expected_lat_range = minimum(getindex.(expected_points, 2)):maximum(getindex.(expected_points, 2))
    expected_mask_sub = expected_mask[expected_lon_range, expected_lat_range]
    expected_arr = Array(view(rotated_cube, expected_lon_range, expected_lat_range)) .* Float32.(expected_mask_sub)

    @test Set(keys(rotated_sub.cubes)) == Set([:pr, :lon, :lat, :rotated_pole])
    @test size(rotated_sub.pr) == size(expected_arr)

    actual_arr = Array(rotated_sub.pr)
    @test all((isnan.(expected_arr) .& isnan.(actual_arr)) .| (expected_arr .== actual_arr))
    @test size(rotated_sub.lon) == size(rotated_sub.pr)
    @test size(rotated_sub.lat) == size(rotated_sub.pr)
    @test rotated_sub.rotated_pole.properties == rp_cube.properties
end
