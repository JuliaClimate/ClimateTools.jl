using ClimateTools
using Test
using YAXArrays
using DimensionalData
using Dates

function _assert_geomakie_fallback(thunk)
    err = try
        thunk()
        nothing
    catch caught
        caught
    end

    @test err isa ArgumentError
    message = sprint(showerror, err)
    @test occursin("GeoMakie", message)
    @test occursin("backend", message)
end

@testset "Plotting wrapper fallback" begin
    if !ClimateTools._ensure_plotting_extension_loaded()
        lon = [0.0, 30.0]
        lat = [40.0, 50.0]
        times = [Date(2020, 1, 1), Date(2020, 1, 2)]
        map_cube = YAXArray(
            (Dim{:longitude}(lon), Dim{:latitude}(lat), Dim{:time}(times)),
            reshape(Float64.(1:8), 2, 2, 2),
        )
        ts_cube = YAXArray((Dim{:time}(times),), Float64[1, 2])

        _assert_geomakie_fallback(() -> geomap(map_cube; dim=:time, index=1))
        _assert_geomakie_fallback(() -> geomapfacet(map_cube; facetdim=:time, ncols=2))
        _assert_geomakie_fallback(() -> timeseriesplot(ts_cube))
        _assert_geomakie_fallback(() -> statsplot(Float64[1, 2, 3]))
        _assert_geomakie_fallback(() -> robustnessmap(YAXArray((Dim{:longitude}(lon), Dim{:latitude}(lat)), reshape([1.0, 2.0, 3.0, 1.0], 2, 2), Dict("flag_values" => [1, 2, 3], "flag_descriptions" => ["Robust signal", "No change or no signal", "Conflicting signal"]))))
    else
        err = try
            ClimateTools._plotting_extension_error(:geomap)
            nothing
        catch caught
            caught
        end

        @test err isa ArgumentError
        message = sprint(showerror, err)
        @test occursin("GeoMakie", message)
        @test occursin("backend", message)
    end
end