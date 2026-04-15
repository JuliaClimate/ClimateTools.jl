# This test file includes parity fixtures and expected values adapted from
# xclim (https://github.com/Ouranosinc/xclim), Copyright 2018-2023 Ouranos Inc.
# and contributors, and is distributed under the Apache License, Version 2.0.
# The ClimateTools.jl version modifies those cases for Julia and YAXArrays.
# See LICENSES/xclim-APACHE-2.0.txt and LICENSE.md for details.

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

    stats_cube = make_member_time_cube([
        1.0 2.0 3.0;
        2.0 4.0 6.0;
        3.0 NaN 9.0;
        4.0 8.0 12.0;
    ])

    stats = ensemble_mean_std_max_min(stats_cube)
    @test vector_values(stats.mean) ≈ [2.5, 14 / 3, 7.5]
    @test vector_values(stats.stdev) ≈ [sqrt(1.25), sqrt(56 / 9 / 3), sqrt(11.25)] atol=1e-10
    @test vector_values(stats.max) == [4.0, 8.0, 12.0]
    @test vector_values(stats.min) == [1.0, 2.0, 3.0]

    weighted_stats = ensemble_mean_std_max_min(stats_cube; weights=[1.0, 0.1, 3.5, 5.0])
    @test vector_values(weighted_stats.mean) ≈ [3.1979166666666665, 7.171171171171171, 9.197916666666666] atol=1e-12
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
end