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
end
