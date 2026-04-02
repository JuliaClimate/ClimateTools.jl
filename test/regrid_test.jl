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
