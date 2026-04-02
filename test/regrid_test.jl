using Test
using ClimateTools
using YAXArrays
using DimensionalData
using Dates

@testset "Regridder rectilinear workflow" begin
    src_lon = collect(0.0:1.0:4.0)
    src_lat = collect(10.0:1.0:13.0)
    dst_lon = collect(0.5:1.0:3.5)
    dst_lat = collect(10.5:1.0:12.5)

    source_field = [2.0 * lon + 3.0 * lat for lon in src_lon, lat in src_lat]
    source2d = YAXArray((Dim{:longitude}(src_lon), Dim{:latitude}(src_lat)), source_field)
    target_grid = YAXArray((Dim{:longitude}(dst_lon), Dim{:latitude}(dst_lat)), zeros(length(dst_lon), length(dst_lat)))

    regridder = Regridder(source2d, target_grid; method="bilinear")
    out2d = regridder(source2d)

    expected2d = [2.0 * lon + 3.0 * lat for lon in dst_lon, lat in dst_lat]
    @test size(out2d) == size(expected2d)
    @test Base.maximum(abs.(Array(out2d) .- expected2d)) < 1e-6

    times = collect(Date(2001, 1, 1):Day(1):Date(2001, 1, 2))
    source3d = YAXArray(
        (Dim{:longitude}(src_lon), Dim{:latitude}(src_lat), Dim{:time}(times)),
        cat(source_field, source_field .+ 1.0; dims=3),
    )

    out3d = regridder(source3d)
    @test size(out3d) == (length(dst_lon), length(dst_lat), 2)
    @test Base.maximum(abs.(Array(out3d)[:, :, 1] .- expected2d)) < 1e-6
    @test Base.maximum(abs.(Array(out3d)[:, :, 2] .- (expected2d .+ 1.0))) < 1e-6

    out_compat = regrid_cube(source3d, target_grid; method="linear")
    @test Array(out_compat) ≈ Array(out3d)
end

@testset "Regridder nearest and skipna" begin
    src_lon = [0.0, 1.0, 2.0]
    src_lat = [0.0, 1.0, 2.0]
    dst_lon = [0.49, 1.51]
    dst_lat = [0.49, 1.51]

    data = [10.0 20.0 30.0; 40.0 50.0 60.0; 70.0 80.0 90.0]
    source = YAXArray((Dim{:longitude}(src_lon), Dim{:latitude}(src_lat)), data)
    dest = YAXArray((Dim{:longitude}(dst_lon), Dim{:latitude}(dst_lat)), zeros(2, 2))

    nearest_regridder = Regridder(source, dest; method="nearest_s2d")
    out_nearest = nearest_regridder(source)

    # nearest in both axes -> picks source indices (1,1), (3,1), (1,3), (3,3)
    @test Array(out_nearest) == [10.0 30.0; 70.0 90.0]

    data_nan = copy(data)
    data_nan[2, 2] = NaN
    source_nan = YAXArray((Dim{:longitude}(src_lon), Dim{:latitude}(src_lat)), data_nan)

    center_dest = YAXArray((Dim{:longitude}([1.5]), Dim{:latitude}([1.5])), zeros(1, 1))
    bilinear = Regridder(source_nan, center_dest; method="bilinear")

    out_no_skip = regrid(source_nan, bilinear; skipna=false)
    @test isnan(Array(out_no_skip)[1, 1])

    out_skip_strict = regrid(source_nan, bilinear; skipna=true, na_thres=0.2)
    @test isnan(Array(out_skip_strict)[1, 1])

    out_skip_relaxed = regrid(source_nan, bilinear; skipna=true, na_thres=0.3)
    @test isfinite(Array(out_skip_relaxed)[1, 1])
end

@testset "Regridder serialization" begin
    src_lon = reverse(collect(0.0:1.0:4.0))
    src_lat = collect(10.0:1.0:13.0)
    dst_lon = collect(0.5:1.0:3.5)
    dst_lat = collect(10.5:1.0:12.5)

    source_field = [2.0 * lon + 3.0 * lat for lon in src_lon, lat in src_lat]
    source = YAXArray((Dim{:longitude}(src_lon), Dim{:latitude}(src_lat)), source_field)
    dest = YAXArray((Dim{:longitude}(dst_lon), Dim{:latitude}(dst_lat)), zeros(length(dst_lon), length(dst_lat)))

    regridder = Regridder(source, dest; method="bilinear")
    original = regrid(source, regridder)

    mktempdir() do tmpdir
        weights_path = joinpath(tmpdir, "regridder.bin")
        saved_path = save_regridder(weights_path, regridder)
        @test saved_path == weights_path

        restored = load_regridder(weights_path)
        @test restored isa Regridder
        @test Array(regrid(source, restored)) ≈ Array(original)

        shifted_source = YAXArray(
            (Dim{:longitude}(src_lon .+ 0.25), Dim{:latitude}(src_lat)),
            source_field,
        )
        @test_throws ErrorException regrid(shifted_source, restored)
    end
end

@testset "Curvilinear interfaces" begin
    lonvec = [0.0, 1.0, 2.0]
    latvec = [10.0, 11.0, 12.0]
    lon2d = repeat(reshape(lonvec, :, 1), 1, length(latvec))
    lat2d = repeat(reshape(latvec, 1, :), length(lonvec), 1)
    data = [lon + lat for lon in lonvec, lat in latvec]

    out_array = regrid_curvilinear_to_regular(lon2d, lat2d, data, lonvec, latvec; method="linear")
    @test size(out_array) == size(data)
    @test out_array ≈ data

    source2d = YAXArray((Dim{:longitude}(lonvec), Dim{:latitude}(latvec)), data)
    target2d = YAXArray((Dim{:longitude}(lonvec), Dim{:latitude}(latvec)), zeros(length(lonvec), length(latvec)))
    out_yax = regrid_curvilinear_to_regular(source2d, target2d; method="linear")
    @test Array(out_yax) ≈ data

    # Rotated path sanity check
    lonrot = [0.0 10.0 20.0; 5.0 15.0 25.0; 10.0 20.0 30.0]
    latrot = [-10.0 0.0 10.0; -5.0 5.0 15.0; 0.0 10.0 20.0]
    rotated_data = [1.0 2.0 3.0; 1.5 2.5 3.5; 2.0 3.0 4.0]

    londest = collect(-180.0:120.0:180.0)
    latdest = collect(-20.0:15.0:25.0)

    out_rot = regrid_rotated_curvilinear_to_regular(
        lonrot,
        latrot,
        rotated_data,
        londest,
        latdest;
        grid_north_longitude=10.0,
        grid_north_latitude=45.0,
        method="idw",
    )

    @test size(out_rot) == (length(londest), length(latdest))
    @test all(isfinite, out_rot)
end

@testset "Regridder Dataset rotated-pole workflow" begin
    # --- simulate a rotated-pole Dataset (like ClimEx) ---
    # Rotated pole: north pole at (0, 42.5) in geographic coords
    north_pole_lon = 0.0
    north_pole_lat = 42.5

    rlon_vals = collect(0.0:1.0:5.0)
    rlat_vals = collect(-2.0:1.0:3.0)
    nrlon = length(rlon_vals)
    nrlat = length(rlat_vals)

    # Compute 2D geographic coordinates from the rotated grid
    rlon2d, rlat2d = ClimateTools.ndgrid(rlon_vals, rlat_vals)
    lon2d, lat2d = ClimateTools.rotated_to_geographic(rlon2d, rlat2d, north_pole_lon, north_pole_lat)

    # Build synthetic source data: a linear function of geographic lon
    source_data = Float32.(lon2d .+ lat2d)

    # Build the cubes for the Dataset
    pr_cube = YAXArray(
        (Dim{:rlon}(rlon_vals), Dim{:rlat}(rlat_vals)),
        source_data,
        Dict{String, Any}("grid_mapping" => "rotated_pole"),
    )
    lon_cube = YAXArray(
        (Dim{:rlon}(rlon_vals), Dim{:rlat}(rlat_vals)),
        Float64.(lon2d),
    )
    lat_cube = YAXArray(
        (Dim{:rlon}(rlon_vals), Dim{:rlat}(rlat_vals)),
        Float64.(lat2d),
    )
    rp_cube = YAXArray(
        (Dim{:maxStrlen64}(1:1),),
        [' '],
        Dict{String, Any}(
            "grid_north_pole_longitude" => north_pole_lon,
            "grid_north_pole_latitude" => north_pole_lat,
        ),
    )

    ds = Dataset(pr=pr_cube, lon=lon_cube, lat=lat_cube, rotated_pole=rp_cube)

    # Destination grid covering a subset of the geographic domain
    lon_range = range(Base.minimum(lon2d) + 0.5, Base.maximum(lon2d) - 0.5, length=4)
    lat_range = range(Base.minimum(lat2d) + 0.5, Base.maximum(lat2d) - 0.5, length=4)
    dst_lon = Float64.(collect(lon_range))
    dst_lat = Float64.(collect(lat_range))
    dest_grid = YAXArray(
        (Dim{:longitude}(dst_lon), Dim{:latitude}(dst_lat)),
        zeros(length(dst_lon), length(dst_lat)),
    )

    # --- Test Regridder from Dataset ---
    regridder = Regridder(ds, :pr, dest_grid)
    @test regridder.weights isa ClimateTools.IDWWeights
    @test regridder.weights.nx_out == length(dst_lon)
    @test regridder.weights.ny_out == length(dst_lat)

    out = regrid(ds[:pr], regridder)
    @test size(out) == (length(dst_lon), length(dst_lat))
    # The interpolated field should not be constant (the original bug)
    out_arr = Array(out)
    @test Base.maximum(out_arr) - Base.minimum(out_arr) > 0.1
    # Output should be finite
    @test all(isfinite, out_arr)

    # --- Test callable form ---
    out2 = regridder(ds[:pr])
    @test Array(out2) ≈ out_arr

    # --- Test regrid_cube(Dataset, ...) convenience ---
    out3 = regrid_cube(ds, :pr, dest_grid)
    @test Array(out3) ≈ out_arr

    dest_ds = Dataset(pr=dest_grid)
    out4 = regrid_cube(ds, :pr, dest_ds)
    @test Array(out4) ≈ out_arr

    regridder_from_dataset_dest = Regridder(ds, :pr, dest_ds)
    @test Array(regrid(ds[:pr], regridder_from_dataset_dest)) ≈ out_arr

    times = [Date(2020, 1, 1), Date(2020, 1, 2), Date(2020, 1, 3)]

    dest_aux = YAXArray(
        (Dim{:longitude}(dst_lon), Dim{:latitude}(dst_lat), Dim{:time}(times)),
        zeros(length(dst_lon), length(dst_lat), length(times)),
    )
    dest_multi = Dataset(aux=dest_aux, tas=dest_aux)
    out5 = regrid_cube(ds, :pr, dest_multi)
    @test Array(out5) ≈ out_arr

    # --- Test that plain YAXArray with grid_mapping is rejected ---
    @test_throws ErrorException Regridder(ds[:pr], dest_grid)

    # --- Test with 3D source (rlon × rlat × time) ---
    data_3d = cat(source_data, source_data .+ 1f0, source_data .+ 2f0; dims=3)
    pr3d = YAXArray(
        (Dim{:rlon}(rlon_vals), Dim{:rlat}(rlat_vals), Dim{:time}(times)),
        data_3d,
        Dict{String, Any}("grid_mapping" => "rotated_pole"),
    )

    ds3d = Dataset(pr=pr3d, lon=lon_cube, lat=lat_cube, rotated_pole=rp_cube)
    regridder3d = Regridder(ds3d, :pr, dest_grid)
    out_3d = regrid(ds3d[:pr], regridder3d)

    @test size(out_3d) == (length(dst_lon), length(dst_lat), 3)
    # Time slices should differ
    @test !isapprox(Array(out_3d)[:, :, 1], Array(out_3d)[:, :, 3])

    # --- Test serialization round-trip ---
    mktempdir() do tmpdir
        path = joinpath(tmpdir, "rotated_regridder.bin")
        save_regridder(path, regridder)
        restored = load_regridder(path)
        @test restored.weights isa ClimateTools.IDWWeights
        out_restored = regrid(ds[:pr], restored)
        @test Array(out_restored) ≈ out_arr
    end

    # --- Test skipna with IDW ---
    source_nan = copy(source_data)
    source_nan[3, 3] = NaN32
    pr_nan = YAXArray(
        (Dim{:rlon}(rlon_vals), Dim{:rlat}(rlat_vals)),
        source_nan,
        Dict{String, Any}("grid_mapping" => "rotated_pole"),
    )
    ds_nan = Dataset(pr=pr_nan, lon=lon_cube, lat=lat_cube, rotated_pole=rp_cube)
    regridder_nan = Regridder(ds_nan, :pr, dest_grid)

    out_no_skip = regrid(ds_nan[:pr], regridder_nan; skipna=false)
    out_skip = regrid(ds_nan[:pr], regridder_nan; skipna=true, na_thres=0.5)
    # With skipna=true, more points should be finite
    @test count(isfinite, Array(out_skip)) >= count(isfinite, Array(out_no_skip))
end

@testset "Regridder Dataset without grid_mapping delegates to rectilinear" begin
    # A dataset with regular geographic coordinates should work via Dataset path
    lon = collect(0.0:1.0:4.0)
    lat = collect(10.0:1.0:13.0)
    data = [2.0 * x + 3.0 * y for x in lon, y in lat]

    cube = YAXArray((Dim{:longitude}(lon), Dim{:latitude}(lat)), data)
    ds = Dataset(pr=cube)

    dst_lon = collect(0.5:1.0:3.5)
    dst_lat = collect(10.5:1.0:12.5)
    dest = YAXArray((Dim{:longitude}(dst_lon), Dim{:latitude}(dst_lat)), zeros(length(dst_lon), length(dst_lat)))

    regridder = Regridder(ds, :pr, dest; method="bilinear")
    @test regridder.weights isa ClimateTools.BilinearWeights

    out = regrid(ds[:pr], regridder)
    expected = [2.0 * x + 3.0 * y for x in dst_lon, y in dst_lat]
    @test Base.maximum(abs.(Array(out) .- expected)) < 1e-6

    dest_ds = Dataset(pr=dest)
    out_ds = regrid_cube(ds, :pr, dest_ds; method="bilinear")
    @test Base.maximum(abs.(Array(out_ds) .- expected)) < 1e-6

    aux = YAXArray((Dim{:longitude}(dst_lon), Dim{:latitude}(dst_lat)), fill(1.0, length(dst_lon), length(dst_lat)))
    dest_multi = Dataset(aux=aux, tas=dest)
    out_multi = regrid_cube(ds, :pr, dest_multi; method="bilinear")
    @test Base.maximum(abs.(Array(out_multi) .- expected)) < 1e-6
end
