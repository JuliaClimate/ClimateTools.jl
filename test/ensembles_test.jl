# This test file includes parity fixtures and expected values adapted from
# xclim (https://github.com/Ouranosinc/xclim), Copyright 2018-2023 Ouranos Inc.
# and contributors, and is distributed under the Apache License, Version 2.0.
# The ClimateTools.jl version modifies those cases for Julia and YAXArrays.
# See LICENSES/xclim-APACHE-2.0.txt and LICENSE.md for details.

using ClimateTools
using Test
using YAXArrays
using DimensionalData
using Dates
using NetCDF
using Statistics

@testset "xclim-style ensembles" begin
    function scalar_value(cube)
        return Array(cube)[]
    end

    function vector_values(cube)
        return vec(Float64.(Array(cube)))
    end

    function make_member_time_cube(data)
        times = collect(DateTime(2000, 1, 1):Year(1):(DateTime(2000, 1, 1) + Year(size(data, 2) - 1)))
        return YAXArray((Dim{:realization}(1:size(data, 1)), Dim{:time}(times)), Float64.(data))
    end

    function write_member_dataset(path, tas_values, pr_values)
        lon = [10.0, 20.0]
        time = [1, 2, 3]
        nccreate(path, "tas", "lon", lon, "time", time, atts=Dict("units" => "K"))
        ncwrite(Float64.(tas_values), path, "tas")
        nccreate(path, "pr", "lon", lon, "time", time, atts=Dict("units" => "mm"))
        ncwrite(Float64.(pr_values), path, "pr")
    end

    function write_member_zarr_dataset(path, tas_values, pr_values)
        dims = (Dim{:lon}([10.0, 20.0]), Dim{:time}([1, 2, 3]))
        ds = Dataset(
            tas=YAXArray(dims, Float64.(tas_values)),
            pr=YAXArray(dims, Float64.(pr_values)),
        )
        savedataset(ds, path=path, driver=:zarr, overwrite=true)
    end

    @testset "open_ensemble_dataset" begin
        mktempdir() do tmpdir
            member1_file = joinpath(tmpdir, "member1.nc")
            member2_file = joinpath(tmpdir, "member2.nc")

            member1_tas = reshape(Float64[1, 2, 3, 4, 5, 6], 2, 3)
            member2_tas = reshape(Float64[11, 12, 13, 14, 15, 16], 2, 3)
            member1_pr = reshape(Float64[0.1, 0.2, 0.3, 0.4, 0.5, 0.6], 2, 3)
            member2_pr = reshape(Float64[1.1, 1.2, 1.3, 1.4, 1.5, 1.6], 2, 3)

            write_member_dataset(member1_file, member1_tas, member1_pr)
            write_member_dataset(member2_file, member2_tas, member2_pr)

            ds = open_ensemble_dataset([member1_file, member2_file]; driver=:netcdf)
            @test Set(keys(ds.cubes)) == Set([:tas, :pr])
            @test collect(lookup(ds.tas, :realization)) == [1, 2]

            axis_symbols = Symbol[name(axis) for axis in ds.tas.axes]
            realization_position = findfirst(==(:realization), axis_symbols)
            @test !isnothing(realization_position)

            member1_selectors = ntuple(i -> i == realization_position ? 1 : Colon(), ndims(ds.tas))
            member2_selectors = ntuple(i -> i == realization_position ? 2 : Colon(), ndims(ds.pr))
            @test Array(ds.tas[member1_selectors...]) == member1_tas
            @test Array(ds.pr[member2_selectors...]) == member2_pr

            named_ds = open_ensemble_dataset(
                [member1_file, member2_file];
                realization_dim="member",
                realization_values=["r1", "r2"],
                driver=:netcdf,
            )
            @test collect(lookup(named_ds.tas, :member)) == ["r1", "r2"]

            @test_throws ErrorException open_ensemble_dataset(String[])
            @test_throws ErrorException open_ensemble_dataset([member1_file, member2_file]; realization_values=[1])
        end

        mktempdir() do tmpdir
            member1_store = joinpath(tmpdir, "member1.zarr")
            member2_store = joinpath(tmpdir, "member2.zarr")

            member1_tas = reshape(Float64[2, 4, 6, 8, 10, 12], 2, 3)
            member2_tas = reshape(Float64[3, 6, 9, 12, 15, 18], 2, 3)
            member1_pr = reshape(Float64[0.2, 0.4, 0.6, 0.8, 1.0, 1.2], 2, 3)
            member2_pr = reshape(Float64[0.3, 0.6, 0.9, 1.2, 1.5, 1.8], 2, 3)

            write_member_zarr_dataset(member1_store, member1_tas, member1_pr)
            write_member_zarr_dataset(member2_store, member2_tas, member2_pr)

            zarr_ds = open_ensemble_dataset(
                [member1_store, member2_store];
                realization_values=["r1", "r2"],
                driver=:zarr,
            )
            @test collect(lookup(zarr_ds.tas, :realization)) == ["r1", "r2"]

            axis_symbols = Symbol[name(axis) for axis in zarr_ds.pr.axes]
            realization_position = findfirst(==(:realization), axis_symbols)
            @test !isnothing(realization_position)

            member1_selectors = ntuple(i -> i == realization_position ? 1 : Colon(), ndims(zarr_ds.tas))
            member2_selectors = ntuple(i -> i == realization_position ? 2 : Colon(), ndims(zarr_ds.pr))
            @test Array(zarr_ds.tas[member1_selectors...]) == member1_tas
            @test Array(zarr_ds.pr[member2_selectors...]) == member2_pr
        end
    end

    stats_cube = make_member_time_cube([
        1.0 2.0 3.0;
        2.0 4.0 6.0;
        3.0 NaN 9.0;
        4.0 8.0 12.0;
    ])

    stats = ensemble_mean_std_max_min(stats_cube)
    @test vector_values(stats.mean) ≈ [2.5, 14 / 3, 7.5]
    @test vector_values(stats.stdev) ≈ [sqrt(1.25), sqrt(56 / 9), sqrt(11.25)] atol=1e-10
    @test vector_values(stats.max) == [4.0, 8.0, 12.0]
    @test vector_values(stats.min) == [1.0, 2.0, 3.0]

    weighted_stats = ensemble_mean_std_max_min(stats_cube; weights=[1.0, 0.1, 3.5, 5.0])
    @test vector_values(weighted_stats.mean) ≈ [
        (1.0 * 1.0 + 2.0 * 0.1 + 3.0 * 3.5 + 4.0 * 5.0) / (1.0 + 0.1 + 3.5 + 5.0),
        (2.0 * 1.0 + 4.0 * 0.1 + 8.0 * 5.0) / (1.0 + 0.1 + 5.0),
        (3.0 * 1.0 + 6.0 * 0.1 + 9.0 * 3.5 + 12.0 * 5.0) / (1.0 + 0.1 + 3.5 + 5.0),
    ] atol=1e-12
    @test vector_values(weighted_stats.max) == [4.0, 8.0, 12.0]
    @test vector_values(weighted_stats.min) == [1.0, 2.0, 3.0]

    masked_stats = ensemble_mean_std_max_min(stats_cube; min_members=4)
    @test isnan(vector_values(masked_stats.mean)[2])
    @test !isnan(vector_values(masked_stats.mean)[1])
    @test !isnan(vector_values(masked_stats.mean)[3])

    percentile_split = ensemble_percentiles(stats_cube)
    percentile_cube = ensemble_percentiles(stats_cube; split=false)
    @test vector_values(percentile_split.p10) ≈ [1.3, 2.4, 3.9] atol=1e-12
    @test vector_values(percentile_split.p50) ≈ [2.5, 4.0, 7.5] atol=1e-12
    @test vector_values(percentile_split.p90) ≈ [3.7, 7.2, 11.1] atol=1e-12
    @test vector_values(percentile_cube[1, :]) ≈ vector_values(percentile_split.p10) atol=1e-12
    @test vector_values(percentile_cube[2, :]) ≈ vector_values(percentile_split.p50) atol=1e-12
    @test vector_values(percentile_cube[3, :]) ≈ vector_values(percentile_split.p90) atol=1e-12

    hazen_series = YAXArray((Dim{:realization}(1:4),), Float64[1, 2, 3, 4])
    hazen = ensemble_percentiles(hazen_series; values=[90], method="hazen", split=false)
    @test scalar_value(hazen) ≈ Statistics.quantile([1.0, 2.0, 3.0, 4.0], 0.9; alpha=0.5, beta=0.5)

    weighted_series = YAXArray((Dim{:realization}(1:7),), Float64[1, 1.9, 2.2, 3, 3.7, 4.1, 5])
    weighted_percentiles = ensemble_percentiles(
        weighted_series;
        values=[20, 40, 60, 80],
        weights=[0.25, 0.05, 0.15, 0.25, 0.15, 0.10, 0.05],
        split=false,
    )
    @test vec(Array(weighted_percentiles)) ≈ [1.554595, 2.463784, 3.0, 3.518378] atol=1e-6

    masked_percentiles = ensemble_percentiles(stats_cube; min_members=4)
    @test isnan(vector_values(masked_percentiles.p50)[2])

    dataset_percentiles = ensemble_percentiles(Dataset(tas=stats_cube, pr=YAXArray(stats_cube.axes, 2 .* Array(stats_cube))))
    @test Set(keys(dataset_percentiles.cubes)) == Set([:tas_p10, :tas_p50, :tas_p90, :pr_p10, :pr_p50, :pr_p90])

    criteria_cube = YAXArray(
        (Dim{:realization}(1:2), Dim{:time}(DateTime(2001, 1, 1):Year(1):DateTime(2002, 1, 1)), Dim{:lat}([10.0, 20.0])),
        Float64[
            1 2 3 4;
            5 6 7 8;
        ] |> x -> reshape(x, 2, 2, 2),
    )
    criteria = make_criteria(criteria_cube)
    @test size(criteria) == (2, 4)
    @test length(collect(lookup(criteria, :criteria))) == 4

    criteria_dataset = Dataset(
        tas=criteria_cube,
        pr=YAXArray((Dim{:realization}(1:2), Dim{:time}(DateTime(2001, 1, 1):Year(1):DateTime(2002, 1, 1))), [1.0 NaN; 2.0 NaN]),
    )
    criteria_ds = make_criteria(criteria_dataset)
    @test size(criteria_ds) == (2, 5)

    kkz_cube = YAXArray((Dim{:realization}(1:4), Dim{:criteria}(["c"])), reshape(Float64[0, 2, 5, 10], 4, 1))
    @test kkz_reduce_ensemble(kkz_cube, 3; standardize=false) == [3, 1, 4]
    @test kkz_reduce_ensemble(YAXArray((Dim{:realization}(1:4), Dim{:time}(1:1)), reshape(Float64[0, 2, 5, 10], 4, 1)), 3; standardize=false) == [3, 1, 4]

    delta = YAXArray((Dim{:realization}(1:6),), Float64[-2, 1, -2, -1, 0, 0])
    fracs = robustness_fractions(delta; test="threshold", abs_thresh=1.5)
    @test scalar_value(fracs.changed) == 2 / 6
    @test scalar_value(fracs.changed_positive) == 0.0
    @test scalar_value(fracs.positive) == 1 / 6
    @test scalar_value(fracs.agree) == 3 / 6

    weighted_delta = YAXArray((Dim{:realization}(1:4),), Float64[-2, 1, -2, -1])
    weighted_fracs = robustness_fractions(weighted_delta; test="threshold", abs_thresh=1.5, weights=[4.0, 3.0, 2.0, 1.0])
    @test scalar_value(weighted_fracs.changed) == 0.6
    @test scalar_value(weighted_fracs.positive) == 0.3
    @test scalar_value(weighted_fracs.changed_positive) == 0.0
    @test scalar_value(weighted_fracs.agree) == 0.7

    empty_ref = YAXArray((Dim{:realization}(1:20), Dim{:time}(DateTime(1900, 1, 1):Day(1):DateTime(1900, 1, 10))), fill(NaN, 20, 10))
    empty_fut = YAXArray((Dim{:realization}(1:20), Dim{:time}(DateTime(1900, 1, 1):Day(1):DateTime(1900, 1, 10))), fill(NaN, 20, 10))
    empty_fracs = robustness_fractions(empty_fut, empty_ref; test="ttest")
    @test scalar_value(empty_fracs.changed) == 0.0
    @test scalar_value(empty_fracs.valid) == 0.0

    changed = YAXArray((Dim{:lat}(1:4),), Float64[0.5, 0.8, 1.0, 1.0])
    agree = YAXArray((Dim{:lat}(1:4),), Float64[1.0, 0.5, 0.5, 1.0])
    categories = robustness_categories(changed, agree)
    @test vec(Array(categories)) == [2, 3, 3, 1]
    @test categories.properties["flag_values"] == [1, 2, 3]
    @test categories.properties["flag_meanings"] == "robust_signal no_change_or_no_signal conflicting_signal"

    ref = YAXArray((Dim{:time}(1:6),), Float64[274, 275, 274.5, 276, 274.3, 273.3])
    fut = YAXArray(
        (Dim{:realization}(1:2), Dim{:time}(1:6)),
        Float64[
            277 277.1 278 278.4 278.1 276.9;
            275 275.8 276 275.2 276.2 275.7;
        ],
    )
    @test scalar_value(robustness_coefficient(fut, ref)) ≈ 0.91972477 atol=1e-6
    @test scalar_value(robustness_coefficient(Dataset(tas=fut), Dataset(tas=ref)).tas) ≈ 0.91972477 atol=1e-6

    @testset "dataset ensemble summaries and validations" begin
        static_cube = YAXArray((stats_cube.axes[2],), Float64[10, 20, 30])
        weights_cube = YAXArray((Dim{:realization}(1:4),), Float64[1.0, 0.1, 3.5, 5.0])
        ds = Dataset(
            tas=stats_cube,
            pr=YAXArray(stats_cube.axes, 2 .* Array(stats_cube)),
            static=static_cube,
        )

        dataset_stats = ensemble_mean_std_max_min(ds; weights=weights_cube)
        @test Set(keys(dataset_stats.cubes)) == Set([
            :tas_mean,
            :tas_stdev,
            :tas_max,
            :tas_min,
            :pr_mean,
            :pr_stdev,
            :pr_max,
            :pr_min,
        ])
        @test vector_values(dataset_stats.tas_mean) ≈ vector_values(weighted_stats.mean) atol=1e-12
        @test vector_values(dataset_stats.pr_mean) ≈ 2 .* vector_values(weighted_stats.mean) atol=1e-12

        percentile_dataset = ensemble_percentiles(ds; values=[12.5, 50], weights=weights_cube, split=false)
        @test Set(keys(percentile_dataset.cubes)) == Set([:tas, :pr])
        @test collect(lookup(percentile_dataset.tas, :percentiles)) == [12.5, 50.0]

        custom_percentiles = ensemble_percentiles(stats_cube; values=[12.5])
        custom_dataset_percentiles = ensemble_percentiles(Dataset(tas=stats_cube); values=[12.5])
        @test Set(keys(custom_percentiles.cubes)) == Set([:p12_5])
        @test Set(keys(custom_dataset_percentiles.cubes)) == Set([:tas_p12_5])

        weighted_cube_percentiles = ensemble_percentiles(stats_cube; values=[20, 50], weights=weights_cube, split=false)
        weighted_vector_percentiles = ensemble_percentiles(stats_cube; values=[20, 50], weights=[1.0, 0.1, 3.5, 5.0], split=false)
        @test Array(weighted_cube_percentiles) ≈ Array(weighted_vector_percentiles) atol=1e-12

        @test_throws ErrorException ensemble_mean_std_max_min(stats_cube; min_members=0)
        @test_throws ErrorException ensemble_mean_std_max_min(Dataset(static=static_cube))
        @test_throws ErrorException ensemble_percentiles(stats_cube; method="not-a-method")
        @test_throws ErrorException ensemble_percentiles(Dataset(static=static_cube))
        @test_throws ErrorException ensemble_percentiles(stats_cube; weights=[1.0, 2.0], split=false)
        @test_throws ErrorException Array(ensemble_percentiles(stats_cube; values=[50], weights=weights_cube, method="hazen", split=false))
    end

    @testset "criteria and kkz branches" begin
        criteria_input = YAXArray((Dim{:realization}(1:3),), Float64[1.0, NaN, 3.0])
        one_dimensional_criteria = make_criteria(criteria_input)
        @test size(one_dimensional_criteria) == (3, 1)
        @test collect(lookup(one_dimensional_criteria, :criteria)) == ["value"]

        cube_a = YAXArray((Dim{:realization}(1:2), Dim{:time}(1:2)), Float64[1 2; 3 4])
        cube_b = YAXArray((Dim{:realization}(2:3), Dim{:time}(1:2)), Float64[5 6; 7 8])
        static_cube = YAXArray((Dim{:time}(1:2),), Float64[1, 2])
        @test_throws ErrorException make_criteria(Dataset(a=cube_a, b=cube_b))
        @test_throws ErrorException make_criteria(Dataset(static=static_cube))

        criteria_cube = YAXArray(
            (Dim{:realization}(1:4), Dim{:criteria}(["a", "b"])),
            Float64[
                0 0;
                1 0;
                0 1;
                3 3;
            ],
        )
        for method in ["cityblock", "chebyshev", "cosine", "seuclidean", "mahalanobis"]
            selected = kkz_reduce_ensemble(criteria_cube, 2; dist_method=method, standardize=false)
            @test length(selected) == 2
            @test length(unique(selected)) == 2
            @test all((1 .<= selected) .& (selected .<= 4))
        end

        kkz_dataset = Dataset(tas=YAXArray((Dim{:realization}(1:4), Dim{:time}(1:1)), reshape(Float64[0, 2, 5, 10], 4, 1)))
        @test kkz_reduce_ensemble(kkz_dataset, 3; standardize=false) == kkz_reduce_ensemble(make_criteria(kkz_dataset), 3; standardize=false)

        nan_criteria_cube = YAXArray((Dim{:realization}(1:2), Dim{:criteria}(["a"])), reshape(Float64[1.0, NaN], 2, 1))
        @test_throws ErrorException kkz_reduce_ensemble(nan_criteria_cube, 1)
        @test_throws ErrorException kkz_reduce_ensemble(criteria_cube, 5)
    end

    @testset "robustness branches" begin
        lats = [10.0, 20.0]
        times = 1:5

        ref_data = zeros(Float64, 2, 2, 5)
        ref_data[1, 1, :] = [1.0, 2.0, 1.0, 2.0, 1.0]
        ref_data[1, 2, :] = [1.5, 2.5, 1.5, 2.5, 1.5]
        ref_data[2, 1, :] = [1.0, 1.5, 1.0, 1.5, 1.0]
        ref_data[2, 2, :] = [0.5, 1.0, 0.5, 1.0, 0.5]

        fut_data = zeros(Float64, 2, 2, 5)
        fut_data[1, 1, :] = [3.0, 4.0, 3.0, 4.0, 3.0]
        fut_data[1, 2, :] = [2.5, 3.5, 2.5, 3.5, 2.5]
        fut_data[2, 1, :] = [1.2, 1.6, 1.1, 1.5, 1.0]
        fut_data[2, 2, :] = [0.8, 1.2, 0.7, 1.1, 0.6]

        future_cube = YAXArray((Dim{:lat}(lats), Dim{:realization}(1:2), Dim{:time}(times)), fut_data)
        reference_cube = YAXArray((Dim{:lat}(lats), Dim{:realization}(1:2), Dim{:time}(times)), ref_data)

        for test_name in ["ttest", "welch-ttest", "mannwhitney-utest", "brownforsythe-test"]
            fractions = robustness_fractions(future_cube, reference_cube; test=test_name)
            @test haskey(fractions.cubes, :pvals)
            @test size(Array(fractions.pvals)) == (2, 2)
            @test all(isfinite, vec(Array(fractions.pvals)))
            @test all((0 .<= vec(Array(fractions.changed))) .& (vec(Array(fractions.changed)) .<= 1))
        end

        invalid_future_data = copy(fut_data)
        invalid_future_data[2, 2, 4:5] .= NaN
        invalid_future_cube = YAXArray((Dim{:lat}(lats), Dim{:realization}(1:2), Dim{:time}(times)), invalid_future_data)
        invalid_fractions = robustness_fractions(invalid_future_cube, reference_cube; test="welch-ttest", invalid=(n=4,))
        @test vec(Array(invalid_fractions.valid)) == [1.0, 0.5]

        delta_cube = YAXArray(
            (Dim{:lat}(lats), Dim{:realization}(1:4)),
            Float64[
                1.0 2.0 3.0 4.0;
                -1.0 0.5 0.0 2.0;
            ],
        )
        rel_threshold_fractions = robustness_fractions(delta_cube; test="threshold", rel_thresh=1.5)
        @test vec(Array(rel_threshold_fractions.changed)) == [0.75, 0.25]
        @test vec(Array(rel_threshold_fractions.agree)) == [1.0, 0.5]

        ipcc_fractions = robustness_fractions(
            future_cube,
            reference_cube;
            test="ipcc-ar6-c",
            ref_pi=Float64.(repeat([0.0, 1.0], 20)),
        )
        @test !haskey(ipcc_fractions.cubes, :pvals)
        @test all(vec(Array(ipcc_fractions.valid)) .== 1.0)

        category_source = Dataset(
            changed=YAXArray((Dim{:lat}(lats),), Float64[0.2, 0.5]),
            agree=YAXArray((Dim{:lat}(lats),), Float64[0.4, 0.7]),
            valid=YAXArray((Dim{:lat}(lats),), Float64[1.0, 0.0]),
        )
        custom_categories = robustness_categories(
            category_source;
            categories=["Exact half", "At most half"],
            ops=[("==", nothing), ("<=", "<=")],
            thresholds=[(0.5, nothing), (0.5, 0.5)],
        )
        @test vec(Array(custom_categories)) == [2, 99]
        @test custom_categories.properties["flag_values"] == [1, 2]
        @test custom_categories.properties["flag_meanings"] == "exact_half at_most_half"

        coefficient = robustness_coefficient(future_cube, reference_cube)
        dataset_coefficient = robustness_coefficient(
            Dataset(tas=future_cube, pr=YAXArray(future_cube.axes, 2 .* Array(future_cube))),
            Dataset(tas=reference_cube, pr=YAXArray(reference_cube.axes, 2 .* Array(reference_cube))),
        )
        @test size(Array(coefficient)) == (2,)
        @test all(isfinite, vec(Array(coefficient)))
        @test Array(dataset_coefficient.tas) ≈ Array(coefficient)
        @test_throws ErrorException robustness_coefficient(Dataset(foo=future_cube), Dataset(bar=reference_cube))
    end
end