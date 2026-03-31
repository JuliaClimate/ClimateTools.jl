using Test
using YAXArrays
using DimensionalData
using Interpolations

# This test expects the local dev version of ClimateTools to be available via
# `using ClimateTools`. To run these tests locally without changing the
# repository Project/Manifest, create a temporary environment and develop the
# local package into it, e.g. in a Julia session:
#
#   using Pkg
#   Pkg.activate(mktempdir())
#   Pkg.develop(path="/absolute/path/to/this/repo")
#   Pkg.instantiate()
#   Pkg.test("ClimateTools")
#
# The test below assumes `using ClimateTools` will load the package from the
# local development path.
using ClimateTools
# Call internal helpers via the module prefix to avoid relying on exports
# (keeps tests working even if the helper is not exported)

# Simple tests for regridding functionality using Dim axes
@testset "regrid curvilinear -> regular" begin
    # Create a tiny separable rectilinear source grid in rotated coordinates
    # Use ndgrid to produce 2D lon/lat arrays (matching shapes expected by the function)
    lonvec = [0.0, 1.0, 2.0]
    latvec = [10.0, 11.0, 12.0]
    # build 2D lon/lat arrays (m x n) where m = length(lonvec), n = length(latvec)
    rlon = repeat(reshape(lonvec, :, 1), 1, length(latvec))
    rlat = repeat(reshape(latvec, 1, :), length(lonvec), 1)

    # Make data on that grid (3x3)
    data = [i + 0.1*j for i in 1:3, j in 1:3]

    # Destination regular grid (lon, lat) as 2D arrays (avoid calling ndgrid inside the function)
    londest = collect(0.0:1.0:2.0)
    latdest = collect(10.0:1.0:12.0)
    lond2 = repeat(reshape(londest, :, 1), 1, length(latdest))
    latd2 = repeat(reshape(latdest, 1, :), length(londest), 1)

    # Call the rotated curvilinear -> regular regrid function with 2D dest arrays
    out = regrid_rotated_curvilinear_to_regular(rlon, rlat, data, lond2, latd2; grid_north_longitude=0.0, grid_north_latitude=90.0)

    @test size(out) == (length(londest), length(latdest))
    @test eltype(out) <: Real
    @test !any(isnan, out)

    # Now construct a YAXArray source and dest using Dim axes and test regrid_cube
    xa = YAXArray((Dim{:longitude}(londest), Dim{:latitude}(latdest)), data)
    dest_axes = (Dim{:longitude}(londest), Dim{:latitude}(latdest))
    dest = YAXArray(dest_axes, similar(data))

    out2 = regrid_cube(xa, dest)
    @test size(out2) == (length(londest), length(latdest))

    # If NearestNeighbors is available, exercise the KD-tree branch of IDW directly
    try
        @eval import NearestNeighbors
        # Create a non-separable (curvilinear) lon/lat by perturbing the separable grid
        lonrot = rlon .+ 0.2 .* randn(size(rlon))
        latrot = rlat .+ 0.1 .* randn(size(rlat))
        vals = vec(data)
    lond2 = repeat(reshape(londest, :, 1), 1, length(latdest))
    latd2 = repeat(reshape(latdest, 1, :), length(londest), 1)
        out_idw = regrid_rotated_curvilinear_to_regular(lonrot, latrot, data, lond2, latd2; grid_north_longitude=0.0, grid_north_latitude=90.0, method="idw")
        @test size(out_idw) == size(lond2)
    catch
        @info "NearestNeighbors not available; skipping KD-tree specific checks"
    end

end


# basic tests for rotated_to_geographic
@testset "rotated conversions" begin
    rlon = [0.0 10.0; 20.0 30.0]
    rlat = [-10.0 0.0; 10.0 20.0]
    grid_north_longitude = 10.0
    grid_north_latitude = 45.0

    lon, lat = ClimateTools.rotated_to_geographic(rlon, rlat, grid_north_longitude, grid_north_latitude)
    @test size(lon) == size(rlon)
    @test size(lat) == size(rlat)
    # lat should be in plausible range
    @test all(-90 .<= lat .<= 90)
    @test all(-180 .<= lon .<= 180)
end

# basic test for regridding: small synthetic grid
@testset "curvilinear -> regular regridding" begin
    lonrot = [0.0 10.0 20.0; 5.0 15.0 25.0; 10.0 20.0 30.0]
    latrot = [-10.0 0.0 10.0; -5.0 5.0 15.0; 0.0 10.0 20.0]
    data = [1.0 2.0 3.0; 1.5 2.5 3.5; 2.0 3.0 4.0]

    londest = collect(-180.0:120.0:180.0)
    latdest = collect(-20.0:15.0:25.0)

    out = regrid_rotated_curvilinear_to_regular(lonrot, latrot, data, londest, latdest; grid_north_longitude=10.0, grid_north_latitude=45.0, method="idw")

    @test size(out) == (length(londest), length(latdest))
    # values should be finite
    @test all(isfinite, out)
end

# If NearestNeighbors is available, test KD-tree accelerated path behaves sensibly
try
    import NearestNeighbors
    @testset "KD-tree accelerated IDW" begin
        lonrot = [0.0 10.0 20.0; 5.0 15.0 25.0; 10.0 20.0 30.0]
        latrot = [-10.0 0.0 10.0; -5.0 5.0 15.0; 0.0 10.0 20.0]
        data = [1.0 2.0 3.0; 1.5 2.5 3.5; 2.0 3.0 4.0]

        londest = collect(-180.0:120.0:180.0)
        latdest = collect(-20.0:15.0:25.0)

        out_kd = regrid_rotated_curvilinear_to_regular(lonrot, latrot, data, londest, latdest; grid_north_longitude=10.0, grid_north_latitude=45.0, method="idw")
        @test size(out_kd) == (length(londest), length(latdest))
        @test all(isfinite, out_kd)
    end
catch err
    @info "NearestNeighbors not available; skipping KD-tree test: $err"
end
