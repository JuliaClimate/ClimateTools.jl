@testset "GeoMakie plotting" begin
    using GeoMakie
    using CairoMakie

    CairoMakie.activate!()

    lon = collect(0.0:30.0:330.0)
    lat = collect(-60.0:30.0:60.0)
    times = [Date(2001, 1, 1), Date(2001, 1, 2), Date(2001, 1, 3)]
    members = 1:4

    data = reshape(Float64.(1:(length(lon) * length(lat) * length(times))), length(lon), length(lat), length(times))
    cube = YAXArray((Dim{:longitude}(lon), Dim{:latitude}(lat), Dim{:time}(times)), data)

    member_data = zeros(Float64, length(lon), length(lat), length(times), length(members))
    for member_index in eachindex(members)
        member_data[:, :, :, member_index] .= data .+ 10 * member_index
    end
    member_cube = YAXArray((Dim{:longitude}(lon), Dim{:latitude}(lat), Dim{:time}(times), Dim{:member}(collect(members))), member_data)

    fig_map = geomap(cube; dim=:time, index=1, colorbar=true, colorbar_label="tas", coastline_width=2.0)
    @test fig_map isa GeoMakie.Makie.Figure
    @test length(fig_map.content) >= 2

    fig_facet = geomapfacet(cube; facetdim=:time, ncols=2, colorbar=true, title="Daily slices")
    @test fig_facet isa GeoMakie.Makie.Figure
    @test length(fig_facet.content) >= 4

    fig_member_facet = geomapfacet(member_cube; facetdim=:member, selectors=(time=1,), ncols=2)
    @test fig_member_facet isa GeoMakie.Makie.Figure

    fig_ts = timeseriesplot(cube; selectors=(longitude=1, latitude=1))
    @test fig_ts isa GeoMakie.Makie.Figure

    fig_ts_reduced = timeseriesplot(cube; reducer=mean)
    @test fig_ts_reduced isa GeoMakie.Makie.Figure

    fig_ts_lines = timeseriesplot(member_cube; selectors=(longitude=1, latitude=1), mode=:lines)
    @test fig_ts_lines isa GeoMakie.Makie.Figure

    fig_ts_ribbon = timeseriesplot(member_cube; selectors=(longitude=1, latitude=1), mode=:mean_ribbon)
    @test fig_ts_ribbon isa GeoMakie.Makie.Figure

    member_stats = ensemble_stats(member_cube; dim="member")
    fig_ts_stats = timeseriesplot(member_stats; selectors=(longitude=1, latitude=1), mode=:stats)
    @test fig_ts_stats isa GeoMakie.Makie.Figure

    fig_stats = statsplot((obs=vec(data[:, :, 1]), corrected=vec(data[:, :, 2])); kind=:hist)
    @test fig_stats isa GeoMakie.Makie.Figure

    fig_box = statsplot(vec(data[:, :, 1]); kind=:boxplot)
    @test fig_box isa GeoMakie.Makie.Figure

    fig_grouped_stats = statsplot(member_cube; groupdim=:member, kind=:boxplot)
    @test fig_grouped_stats isa GeoMakie.Makie.Figure
end

@testset "GeoMakie rotated-pole plotting" begin
    using GeoMakie
    using CairoMakie

    CairoMakie.activate!()

    north_pole_lon = 0.0
    north_pole_lat = 42.5

    rlon_vals = collect(0.0:1.0:5.0)
    rlat_vals = collect(-2.0:1.0:3.0)
    times = [Date(2001, 1, 1), Date(2001, 1, 2)]
    rlon2d, rlat2d = ClimateTools.ndgrid(rlon_vals, rlat_vals)
    lon2d, lat2d = ClimateTools.rotated_to_geographic(rlon2d, rlat2d, north_pole_lon, north_pole_lat)

    source_slice = Float32.(lon2d .+ lat2d)
    source_data = cat(source_slice, source_slice .+ 3; dims=3)
    pr_cube = YAXArray(
        (Dim{:rlon}(rlon_vals), Dim{:rlat}(rlat_vals), Dim{:time}(times)),
        source_data,
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

    ds = Dataset(pr=pr_cube, lon=lon_cube, lat=lat_cube, rotated_pole=rp_cube)

    fig_rot = geomap(ds, :pr; dim=:time, index=1, colorbar=true)
    @test fig_rot isa GeoMakie.Makie.Figure

    fig_rot_facet = geomapfacet(ds, :pr; facetdim=:time, ncols=2)
    @test fig_rot_facet isa GeoMakie.Makie.Figure

    @test_throws ErrorException geomap(pr_cube)
end
