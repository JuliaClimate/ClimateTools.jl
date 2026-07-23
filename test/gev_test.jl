@testset "Extreme value model fitting on cubes" begin
    using Extremes
    using DataFrames
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
    @test size(returnlevel_cube(gev_models), 3) == 7
    @test_throws ArgumentError returnlevel_cube(gev_models; rlevels=[1, 10])
    @test_throws ArgumentError gevfit_cube(gev_cube; locationcov=ones(ntime))

    @testset "Non-stationary GEV covariates" begin
        fremantle = Extremes.dataset("fremantle")
        block_axis = collect(fremantle.Year)
        nblocks = length(block_axis)
        sites = [1, 2]

        response = Matrix{Union{Missing, Float64}}(
            hcat(Float64.(fremantle.SeaLevel), Float64.(fremantle.SeaLevel) .+ 0.05)
        )
        response[5, 1] = missing
        block_cube = YAXArray(
            (Dim{:site}(sites), Dim{:time}(block_axis)),
            permutedims(response),
        )

        year_values = Matrix{Union{Missing, Float64}}(
            hcat(Float64.(block_axis), Float64.(block_axis) .- 1900)
        )
        year_values[7, 1] = missing
        year_cube = YAXArray(
            (Dim{:time}(block_axis), Dim{:site}(sites)),
            year_values,
        )
        soi = Float64.(fremantle.SOI)

        fit = gevfit_cube(
            block_cube;
            covariates=(year=year_cube, soi=soi),
            locationcovid=[:year, :soi],
            logscalecovid=[:year, :soi],
        )

        @test fit.fit_kind == :gev
        @test Tuple(name.(fit.models.axes)) == (:site,)
        @test name(fit.fit_axis) == :time
        @test Tuple(name.(fit.valid_indices.axes)) == (:site,)

        valid_site1 = setdiff(collect(1:nblocks), [5, 7])
        direct_site1 = Extremes.gevfit(
            DataFrame(
                value=Float64.(response[valid_site1, 1]),
                year=Float64.(year_values[valid_site1, 1]),
                soi=soi[valid_site1],
            ),
            :value;
            locationcovid=[:year, :soi],
            logscalecovid=[:year, :soi],
        )
        direct_site2 = Extremes.gevfit(
            DataFrame(
                value=Float64.(response[:, 2]),
                year=Float64.(year_values[:, 2]),
                soi=soi,
            ),
            :value;
            locationcovid=[:year, :soi],
            logscalecovid=[:year, :soi],
        )

        model_array = Array(fit.models)
        @test model_array[1].θ̂ ≈ direct_site1.θ̂
        @test model_array[2].θ̂ ≈ direct_site2.θ̂
        @test Array(fit.valid_indices)[1] == valid_site1
        @test Array(fit.valid_indices)[2] == collect(1:nblocks)

        effective_levels = returnlevel_cube(fit; rlevels=[10, 100])
        effective_array = Array(effective_levels)
        @test Tuple(name.(effective_levels.axes)) == (:site, :time, :rlevels)

        invalid_site1 = setdiff(collect(1:nblocks), valid_site1)
        for (level_index, return_period) in pairs([10, 100])
            expected = Extremes.returnlevel(direct_site1, return_period).value
            actual = effective_array[1, :, level_index]
            @test all(ismissing, actual[invalid_site1])
            @test Float64[value for value in actual if !ismissing(value)] ≈ expected
        end

        shape_fit = gevfit_cube(
            block_cube;
            covariates=(soi=soi,),
            shapecovid=[:soi],
        )
        direct_shape = Extremes.gevfit(
            DataFrame(value=Float64.(response[:, 2]), soi=soi),
            :value;
            shapecovid=[:soi],
        )
        @test Array(shape_fit.models)[2].θ̂ ≈ direct_shape.θ̂

        invalid_fit = gevfit_cube(
            gev_cube;
            covariates=(trend=collect(1.0:ntime),),
            locationcovid=[:trend],
        )
        @test ismissing(Array(invalid_fit.models)[2, 2])
        @test isempty(Array(invalid_fit.valid_indices)[2, 2])
        @test all(ismissing, Array(returnlevel_cube(invalid_fit; rlevels=[2, 5]))[2, 2, :, :])

        save_path = tempname()
        open(save_path, "w") do io
            serialize(io, fit)
        end
        reloaded_fit = open(save_path, "r") do io
            deserialize(io)
        end
        @test isequal(
            Array(returnlevel_cube(reloaded_fit; rlevels=[10, 100])),
            effective_array,
        )

        missing_dim_covariate = YAXArray((Dim{:site}(sites),), Float64.(sites))
        extra_dim_covariate = YAXArray(
            (Dim{:time}(block_axis), Dim{:member}([1])),
            reshape(Float64.(block_axis), nblocks, 1),
        )
        mismatched_time_covariate = YAXArray(
            (Dim{:time}(block_axis .+ 1),),
            Float64.(block_axis),
        )

        @test_throws ArgumentError gevfit_cube(block_cube; covariates=(year=Float64.(block_axis),), locationcovid=[:unknown])
        @test_throws ArgumentError gevfit_cube(block_cube; covariates=(year=Float64.(block_axis), soi=soi), locationcovid=[:year])
        @test_throws ArgumentError gevfit_cube(block_cube; covariates=(year=Float64.(block_axis[1:end-1]),), locationcovid=[:year])
        @test_throws ArgumentError gevfit_cube(block_cube; covariates=(year=missing_dim_covariate,), locationcovid=[:year])
        @test_throws ArgumentError gevfit_cube(block_cube; covariates=(year=extra_dim_covariate,), locationcovid=[:year])
        @test_throws ArgumentError gevfit_cube(block_cube; covariates=(year=mismatched_time_covariate,), locationcovid=[:year])
        @test_throws ArgumentError gevfit_cube(block_cube; covariates=(year=fill("invalid", nblocks),), locationcovid=[:year])
        @test_throws ArgumentError gevfit_cube(block_cube; covariates=(year=Float64.(block_axis),), locationcovid=[:year, :year])
        @test_throws ArgumentError gevfit_cube(block_cube; covariates=(__climatetools_extreme_value__=Float64.(block_axis),), locationcovid=[:__climatetools_extreme_value__])
        @test_throws ArgumentError gevfit_cube(block_cube; dim=(:site, :time), covariates=(year=Float64.(block_axis),), locationcovid=[:year])
        @test_throws ArgumentError gevfit_cube(block_cube; covariates=(year=Float64.(block_axis),))

        underidentified_cube = YAXArray(
            (Dim{:site}([1]), Dim{:time}(1:4)),
            reshape([1.0, 1.2, 1.1, 1.3], 1, :),
        )
        underidentified_fit = gevfit_cube(
            underidentified_cube;
            covariates=(trend=collect(1.0:4.0),),
            locationcovid=[:trend],
            logscalecovid=[:trend],
        )
        @test ismissing(only(Array(underidentified_fit.models)))
        @test_throws ArgumentError returnlevel_cube(fit; rlevels=[1])
    end

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
    @test_throws ArgumentError gpfit_cube(gp_cube; threshold=1.5, logscalecov=ones(ntime))

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
    @test size(returnlevel_cube(gp_with_thresholds), 3) == 8

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