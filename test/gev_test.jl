@testset "Extreme value model fitting on cubes" begin
    using Extremes
    using Serialization
    using StableRNGs

    rng = StableRNG(20260408)
    lon = [10.0, 20.0]
    lat = [45.0, 55.0]
    time = collect(Date(2001, 1, 1):Day(1):Date(2001, 4, 30))
    ntime = length(time)

    gev_data = Array{Float64}(undef, length(lon), length(lat), ntime)
    for i in eachindex(lon), j in eachindex(lat)
        gev_data[i, j, :] = rand(rng, Extremes.GeneralizedExtremeValue(1.5 + 0.25 * i + 0.1 * j, 0.8, 0.1), ntime)
    end
    gev_data[2, 2, :] .= NaN

    gev_cube = YAXArray((Dim{:longitude}(lon), Dim{:latitude}(lat), Dim{:time}(time)), gev_data)
    gev_models = gevfit_cube(gev_cube)
    gev_array = Array(gev_models)

    @test_throws ArgumentError gevfit_cube(gev_cube; minimalvalue=1.0)

    @test Tuple(name.(gev_models.axes)) == (:longitude, :latitude)
    @test gev_array[1, 1] isa Extremes.MaximumLikelihoodAbstractExtremeValueModel
    @test ismissing(gev_array[2, 2])

    direct_gev = Extremes.gevfit(vec(gev_data[1, 1, :]))
    @test gev_array[1, 1].θ̂ ≈ direct_gev.θ̂

    gev_rlevels = returnlevel_cube(gev_models; rlevels=[2, 5, 10])
    gev_rlevel_array = Array(gev_rlevels)
    direct_gev_rlevels = [Extremes.returnlevel(gev_array[1, 1], period).value[1] for period in [2, 5, 10]]

    @test Tuple(name.(gev_rlevels.axes)) == (:longitude, :latitude, :rlevels)
    @test collect(gev_rlevel_array[1, 1, :]) ≈ direct_gev_rlevels

    gp_data = Array{Float64}(undef, length(lon), length(lat), ntime)
    for i in eachindex(lon), j in eachindex(lat)
        gp_data[i, j, :] = 1.5 .+ rand(rng, Extremes.GeneralizedPareto(0.0, 0.8 + 0.1 * i, 0.1), ntime)
    end
    gp_data[2, 2, :] .= NaN

    gp_cube = YAXArray((Dim{:longitude}(lon), Dim{:latitude}(lat), Dim{:time}(time)), gp_data)
    gp_models = gpfit_cube(gp_cube; threshold=1.5)
    gp_array = Array(gp_models)

    @test Tuple(name.(gp_models.axes)) == (:longitude, :latitude)
    @test gp_array[1, 1] isa Extremes.MaximumLikelihoodAbstractExtremeValueModel
    @test ismissing(gp_array[2, 2])

    direct_gp = Extremes.gpfit(vec(gp_data[1, 1, :]) .- 1.5)
    @test gp_array[1, 1].θ̂ ≈ direct_gp.θ̂

    gp_with_thresholds = gpfit_cube(gp_cube; threshold_quantile=0.8, minimalvalue=0.0, return_thresholds=true)
    gp_auto_array = Array(gp_with_thresholds.models)
    threshold_array = Array(gp_with_thresholds.thresholds)

    @test gp_auto_array[1, 2] isa Extremes.MaximumLikelihoodAbstractExtremeValueModel
    @test threshold_array[1, 2] ≈ quantile(vec(gp_data[1, 2, :]), 0.8)
    @test ismissing(threshold_array[2, 2])
    @test gp_with_thresholds.n_observations == ntime
    @test gp_with_thresholds.n_obs_per_block == 365

    gp_rlevels = returnlevel_cube(gp_with_thresholds; rlevels=[2, 5, 10])
    gp_rlevel_array = Array(gp_rlevels)
    direct_gp_rlevels = [
        Extremes.returnlevel(
            gp_auto_array[1, 2],
            threshold_array[1, 2],
            gp_with_thresholds.n_observations,
            gp_with_thresholds.n_obs_per_block,
            period,
        ).value[1] for period in [2, 5, 10]
    ]

    @test Tuple(name.(gp_rlevels.axes)) == (:longitude, :latitude, :rlevels)
    @test collect(gp_rlevel_array[1, 2, :]) ≈ direct_gp_rlevels

    save_path = tempname()
    open(save_path, "w") do io
        serialize(io, gp_with_thresholds)
    end
    gp_reloaded = open(save_path, "r") do io
        deserialize(io)
    end
    gp_rlevels_reloaded = returnlevel_cube(gp_reloaded; rlevels=[2, 5, 10])

    @test isequal(Array(gp_rlevels_reloaded), gp_rlevel_array)
end