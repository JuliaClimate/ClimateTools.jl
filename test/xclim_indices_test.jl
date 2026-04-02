@testset "xclim-style indices" begin
    function scalar_timeseries(cube)
        return vec(Array(cube))
    end

    function make_cube(values, dates)
        return YAXArray((Dim{:longitude}(1:1), Dim{:latitude}(1:1), Dim{:time}(dates)), reshape(values, 1, 1, :))
    end

    function grouped(values, dates, freq)
        if freq == "YS"
            keys = year.(dates)
        else
            keys = yearmonth.(dates)
        end
        return [values[findall(==(key), keys)] for key in unique(keys)]
    end

    function rolling_sum(values, window)
        out = fill(NaN, length(values))
        for i in eachindex(values)
            first_index = i - window + 1
            if first_index >= 1
                out[i] = sum(values[first_index:i])
            end
        end
        return out
    end

    function count_matches(values, comparator, threshold)
        return count(value -> comparator(value, threshold), values)
    end

    function count_pair_matches(first_values, second_values, predicate)
        total = 0
        for (first_value, second_value) in zip(first_values, second_values)
            if !(isnan(first_value) || isnan(second_value)) && predicate(first_value, second_value)
                total += 1
            end
        end
        return total
    end

    function count_matches_against(values, thresholds, comparator)
        total = 0
        for (value, threshold) in zip(values, thresholds)
            if !(isnan(value) || isnan(threshold)) && comparator(value, threshold)
                total += 1
            end
        end
        return total
    end

    function longest_run(values, predicate)
        best = 0
        current = 0
        for value in values
            if predicate(value)
                current += 1
                best = max(best, current)
            else
                current = 0
            end
        end
        return best
    end

    function threshold_sum(values, threshold, comparator, difference)
        total = 0.0
        has_valid = false

        for value in values
            if isnan(value)
                continue
            end

            has_valid = true
            if comparator(value, threshold)
                total += difference(value, threshold)
            end
        end

        return has_valid ? total : NaN
    end

    function run_length_stats(values, predicate, window)
        event_count = 0
        total_days = 0
        max_length = 0
        current = 0
        has_valid = false

        for value in values
            if isnan(value)
                if current >= window
                    event_count += 1
                    total_days += current
                    max_length = max(max_length, current)
                end
                current = 0
                continue
            end

            has_valid = true
            if predicate(value)
                current += 1
            else
                if current >= window
                    event_count += 1
                    total_days += current
                    max_length = max(max_length, current)
                end
                current = 0
            end
        end

        if current >= window
            event_count += 1
            total_days += current
            max_length = max(max_length, current)
        end

        return has_valid ? (event_count, total_days, max_length) : (NaN, NaN, NaN)
    end

    function pair_run_length_stats(first_values, second_values, predicate, window)
        event_count = 0
        total_days = 0
        max_length = 0
        current = 0
        has_valid = false

        for (first_value, second_value) in zip(first_values, second_values)
            if isnan(first_value) || isnan(second_value)
                if current >= window
                    event_count += 1
                    total_days += current
                    max_length = max(max_length, current)
                end
                current = 0
                continue
            end

            has_valid = true
            if predicate(first_value, second_value)
                current += 1
            else
                if current >= window
                    event_count += 1
                    total_days += current
                    max_length = max(max_length, current)
                end
                current = 0
            end
        end

        if current >= window
            event_count += 1
            total_days += current
            max_length = max(max_length, current)
        end

        return has_valid ? (event_count, total_days, max_length) : (NaN, NaN, NaN)
    end

    function spell_total_length_against(values, thresholds, comparator, window)
        total = 0
        current = 0
        has_valid = false

        for (value, threshold) in zip(values, thresholds)
            if isnan(value) || isnan(threshold)
                if current >= window
                    total += current
                end
                current = 0
                continue
            end

            has_valid = true
            if comparator(value, threshold)
                current += 1
            else
                if current >= window
                    total += current
                end
                current = 0
            end
        end

        if current >= window
            total += current
        end

        return has_valid ? total : NaN
    end

    dates = collect(DateTime(2000, 1, 1):Day(1):DateTime(2001, 12, 31))
    tasmax_pattern_2000 = [10.0, 12.0, 14.0, 16.0, 18.0]
    tasmax_pattern_2001 = [20.0, 22.0, 24.0, 26.0, 28.0]
    dtr_pattern_2000 = [5.0, 7.0, 9.0, 11.0, 13.0]
    dtr_pattern_2001 = [6.0, 8.0, 10.0, 12.0, 14.0]
    pr_pattern_2000 = [0.0, 0.5, 1.0, 2.0, 5.0]
    pr_pattern_2001 = [0.0, 1.0, 3.0, 4.0, 6.0]

    tasmax_values = Float64[]
    tasmin_values = Float64[]
    pr_values = Float64[]
    for date in dates
        idx = mod(dayofyear(date) - 1, 5) + 1
        if year(date) == 2000
            push!(tasmax_values, tasmax_pattern_2000[idx])
            push!(tasmin_values, tasmax_pattern_2000[idx] - dtr_pattern_2000[idx])
            push!(pr_values, pr_pattern_2000[idx])
        else
            push!(tasmax_values, tasmax_pattern_2001[idx])
            push!(tasmin_values, tasmax_pattern_2001[idx] - dtr_pattern_2001[idx])
            push!(pr_values, pr_pattern_2001[idx])
        end
    end

    tasmax = make_cube(tasmax_values, dates)
    tasmin = make_cube(tasmin_values, dates)
    pr = make_cube(pr_values, dates)
    dtr_values = tasmax_values .- tasmin_values

    expected_tx_max = [maximum(group) for group in grouped(tasmax_values, dates, "YS")]
    expected_tx_min = [minimum(group) for group in grouped(tasmax_values, dates, "YS")]
    expected_tn_max = [maximum(group) for group in grouped(tasmin_values, dates, "YS")]
    expected_tn_min = [minimum(group) for group in grouped(tasmin_values, dates, "YS")]
    expected_dtr_mean = [mean(group) for group in grouped(dtr_values, dates, "YS")]
    expected_dtr_max = [maximum(group) for group in grouped(dtr_values, dates, "YS")]
    expected_dtr_variability = [mean(abs.(Base.diff(group))) for group in grouped(dtr_values, dates, "YS")]
    expected_rx1day = [maximum(group) for group in grouped(pr_values, dates, "YS")]
    rolling3 = rolling_sum(pr_values, 3)
    expected_rx3day = [maximum(filter(!isnan, group)) for group in grouped(rolling3, dates, "YS")]
    expected_sdii = [mean(filter(x -> x >= 1.0, group)) for group in grouped(pr_values, dates, "YS")]

    @test scalar_timeseries(tx_max(tasmax)) == expected_tx_max
    @test scalar_timeseries(tx_min(tasmax)) == expected_tx_min
    @test scalar_timeseries(tn_max(tasmin)) == expected_tn_max
    @test scalar_timeseries(tn_min(tasmin)) == expected_tn_min
    @test scalar_timeseries(daily_temperature_range(tasmin, tasmax)) ≈ expected_dtr_mean
    @test scalar_timeseries(daily_temperature_range(tasmin, tasmax; op="max")) == expected_dtr_max
    @test scalar_timeseries(daily_temperature_range_variability(tasmin, tasmax)) ≈ expected_dtr_variability
    @test scalar_timeseries(extreme_temperature_range(tasmin, tasmax)) == expected_tx_max .- expected_tn_min
    @test scalar_timeseries(max_1day_precipitation_amount(pr)) == expected_rx1day
    @test scalar_timeseries(max_n_day_precipitation_amount(pr; window=3)) == expected_rx3day
    @test scalar_timeseries(daily_pr_intensity(pr)) ≈ expected_sdii

    short_dates = collect(DateTime(2000, 1, 1):Day(1):DateTime(2000, 2, 4))
    short_values = collect(1.0:length(short_dates))
    short_cube = make_cube(short_values, short_dates)

    @test scalar_timeseries(tx_max(short_cube; freq="MS")) == [31.0, 35.0]
    @test scalar_timeseries(max_n_day_precipitation_amount(short_cube; window=3, freq="MS")) == [90.0, 102.0]
    @test_throws ErrorException tx_max(short_cube; freq="QS")

    tas_values = (tasmax_values .+ tasmin_values) ./ 2
    tas = make_cube(tas_values, dates)
    expected_tg_max = [maximum(group) for group in grouped(tas_values, dates, "YS")]
    expected_tg_mean = [mean(group) for group in grouped(tas_values, dates, "YS")]
    expected_tg_min = [minimum(group) for group in grouped(tas_values, dates, "YS")]

    @test scalar_timeseries(tg_max(tas)) == expected_tg_max
    @test scalar_timeseries(tg_mean(tas)) ≈ expected_tg_mean
    @test scalar_timeseries(tg_min(tas)) == expected_tg_min
    @test scalar_timeseries(tg_mean(short_cube; freq="MS")) == [16.0, 33.5]

    count_tasmax_pattern_2000 = [24.0, 25.0, 26.0, 30.0, 31.0]
    count_tasmax_pattern_2001 = [20.0, 29.0, 30.0, 31.0, 32.0]
    count_tasmin_pattern_2000 = [19.0, 20.0, 21.0, 22.0, 23.0]
    count_tasmin_pattern_2001 = [18.0, 21.0, 22.0, 23.0, 24.0]
    ice_tasmax_pattern_2000 = [-2.0, -1.0, 0.0, 1.0, 2.0]
    ice_tasmax_pattern_2001 = [-5.0, -1.0, 1.0, 3.0, 5.0]

    count_tasmax_values = Float64[]
    count_tasmin_values = Float64[]
    ice_tasmax_values = Float64[]
    for date in dates
        idx = mod(dayofyear(date) - 1, 5) + 1
        if year(date) == 2000
            push!(count_tasmax_values, count_tasmax_pattern_2000[idx])
            push!(count_tasmin_values, count_tasmin_pattern_2000[idx])
            push!(ice_tasmax_values, ice_tasmax_pattern_2000[idx])
        else
            push!(count_tasmax_values, count_tasmax_pattern_2001[idx])
            push!(count_tasmin_values, count_tasmin_pattern_2001[idx])
            push!(ice_tasmax_values, ice_tasmax_pattern_2001[idx])
        end
    end

    count_tasmax = make_cube(count_tasmax_values, dates)
    count_tasmin = make_cube(count_tasmin_values, dates)
    ice_tasmax = make_cube(ice_tasmax_values, dates)

    expected_hot_days = [count_matches(group, >, 25.0) for group in grouped(count_tasmax_values, dates, "YS")]
    expected_warm_days = [count_matches(group, >, 30.0) for group in grouped(count_tasmax_values, dates, "YS")]
    expected_warm_nights = [count_matches(group, >, 22.0) for group in grouped(count_tasmin_values, dates, "YS")]
    expected_ice_days = [count_matches(group, <, 0.0) for group in grouped(ice_tasmax_values, dates, "YS")]

    @test scalar_timeseries(hot_days(count_tasmax)) == expected_hot_days
    @test scalar_timeseries(warm_day_frequency(count_tasmax)) == expected_warm_days
    @test scalar_timeseries(warm_night_frequency(count_tasmin)) == expected_warm_nights
    @test scalar_timeseries(ice_days(ice_tasmax)) == expected_ice_days

    expected_dry_days = [count_matches(group, <, 0.2) for group in grouped(pr_values, dates, "YS")]
    expected_wet_days = [count_matches(group, >=, 1.0) for group in grouped(pr_values, dates, "YS")]
    expected_wet_prop = expected_wet_days ./ [length(group) for group in grouped(pr_values, dates, "YS")]
    expected_mcdd = [longest_run(group, value -> value < 1.0) for group in grouped(pr_values, dates, "YS")]
    expected_mcwd = [longest_run(group, value -> value > 1.0) for group in grouped(pr_values, dates, "YS")]

    @test scalar_timeseries(dry_days(pr)) == expected_dry_days
    @test scalar_timeseries(wetdays(pr)) == expected_wet_days
    @test scalar_timeseries(wetdays_prop(pr)) ≈ expected_wet_prop
    @test scalar_timeseries(maximum_consecutive_dry_days(pr)) == expected_mcdd
    @test scalar_timeseries(maximum_consecutive_wet_days(pr)) == expected_mcwd
    @test scalar_timeseries(wetdays(short_cube; thresh=30.0, freq="MS")) == [2.0, 4.0]

    threshold_tas_pattern_2000 = [-2.0, 0.0, 4.0, 10.0, 12.0]
    threshold_tas_pattern_2001 = [-1.0, 5.0, 7.0, 11.0, 14.0]
    frost_tasmin_pattern_2000 = [-3.0, -1.0, 0.0, 2.0, 4.0]
    frost_tasmin_pattern_2001 = [-5.0, 0.0, 1.0, 3.0, 5.0]

    threshold_tas_values = Float64[]
    frost_tasmin_values = Float64[]
    for date in dates
        idx = mod(dayofyear(date) - 1, 5) + 1
        if year(date) == 2000
            push!(threshold_tas_values, threshold_tas_pattern_2000[idx])
            push!(frost_tasmin_values, frost_tasmin_pattern_2000[idx])
        else
            push!(threshold_tas_values, threshold_tas_pattern_2001[idx])
            push!(frost_tasmin_values, frost_tasmin_pattern_2001[idx])
        end
    end

    threshold_tas = make_cube(threshold_tas_values, dates)
    frost_tasmin = make_cube(frost_tasmin_values, dates)

    expected_frost_days = [count_matches(group, <, 0.0) for group in grouped(frost_tasmin_values, dates, "YS")]
    expected_tg_above = [count_matches(group, >, 10.0) for group in grouped(threshold_tas_values, dates, "YS")]
    expected_tg_below = [count_matches(group, <, 10.0) for group in grouped(threshold_tas_values, dates, "YS")]
    expected_tn_above = [count_matches(group, >, 20.0) for group in grouped(count_tasmin_values, dates, "YS")]
    expected_tn_below = [count_matches(group, <, 20.0) for group in grouped(count_tasmin_values, dates, "YS")]
    expected_tx_above = [count_matches(group, >, 25.0) for group in grouped(count_tasmax_values, dates, "YS")]
    expected_tx_below = [count_matches(group, <, 25.0) for group in grouped(count_tasmax_values, dates, "YS")]
    expected_gdd = [threshold_sum(group, 4.0, >, (value, threshold) -> value - threshold) for group in grouped(threshold_tas_values, dates, "YS")]
    expected_hdd = [threshold_sum(group, 17.0, <, (value, threshold) -> threshold - value) for group in grouped(threshold_tas_values, dates, "YS")]
    expected_cdd = [threshold_sum(group, 10.0, >, (value, threshold) -> value - threshold) for group in grouped(threshold_tas_values, dates, "YS")]
    expected_mcfd = [longest_run(group, value -> value < 0.0) for group in grouped(frost_tasmin_values, dates, "YS")]
    expected_mcffd = [longest_run(group, value -> value >= 0.0) for group in grouped(frost_tasmin_values, dates, "YS")]
    expected_mctx = [longest_run(group, value -> value > 25.0) for group in grouped(count_tasmax_values, dates, "YS")]

    @test scalar_timeseries(frost_days(frost_tasmin)) == expected_frost_days
    @test scalar_timeseries(tg_days_above(threshold_tas)) == expected_tg_above
    @test scalar_timeseries(tg_days_below(threshold_tas)) == expected_tg_below
    @test scalar_timeseries(tn_days_above(count_tasmin)) == expected_tn_above
    @test scalar_timeseries(tn_days_below(count_tasmin; thresh=20.0)) == expected_tn_below
    @test scalar_timeseries(tx_days_above(count_tasmax)) == expected_tx_above
    @test scalar_timeseries(tx_days_below(count_tasmax)) == expected_tx_below
    @test scalar_timeseries(growing_degree_days(threshold_tas)) ≈ expected_gdd
    @test scalar_timeseries(heating_degree_days(threshold_tas)) ≈ expected_hdd
    @test scalar_timeseries(cooling_degree_days(threshold_tas; thresh=10.0)) ≈ expected_cdd
    @test scalar_timeseries(maximum_consecutive_frost_days(frost_tasmin)) == expected_mcfd
    @test scalar_timeseries(maximum_consecutive_frost_free_days(frost_tasmin)) == expected_mcffd
    @test scalar_timeseries(maximum_consecutive_tx_days(count_tasmax)) == expected_mctx

    cold_spell_pattern_2000 = [-12.0, -11.0, -9.0, -12.0, -13.0, -14.0, -8.0, -11.0, -12.0, -9.0]
    cold_spell_pattern_2001 = [-15.0, -12.0, -11.0, -5.0, -4.0, -13.0, -14.0, -15.0, -16.0, -1.0]
    hot_spell_pattern_2000 = [28.0, 31.0, 32.0, 29.0, 33.0, 34.0, 35.0, 20.0, 31.0, 32.0]
    hot_spell_pattern_2001 = [31.0, 32.0, 33.0, 20.0, 29.0, 34.0, 35.0, 36.0, 37.0, 10.0]

    cold_spell_values = Float64[]
    hot_spell_values = Float64[]
    for date in dates
        idx = mod(dayofyear(date) - 1, 10) + 1
        if year(date) == 2000
            push!(cold_spell_values, cold_spell_pattern_2000[idx])
            push!(hot_spell_values, hot_spell_pattern_2000[idx])
        else
            push!(cold_spell_values, cold_spell_pattern_2001[idx])
            push!(hot_spell_values, hot_spell_pattern_2001[idx])
        end
    end

    cold_spell_cube = make_cube(cold_spell_values, dates)
    hot_spell_cube = make_cube(hot_spell_values, dates)

    expected_cold_spell_stats = [run_length_stats(group, value -> value < -10.0, 3) for group in grouped(cold_spell_values, dates, "YS")]
    expected_hot_spell_stats = [run_length_stats(group, value -> value > 30.0, 3) for group in grouped(hot_spell_values, dates, "YS")]
    expected_heat_wave_index = [run_length_stats(group, value -> value > 25.0, 5)[2] for group in grouped(hot_spell_values, dates, "YS")]

    @test scalar_timeseries(cold_spell_days(cold_spell_cube; window=3)) == [stats[2] for stats in expected_cold_spell_stats]
    @test scalar_timeseries(cold_spell_frequency(cold_spell_cube; window=3)) == [stats[1] for stats in expected_cold_spell_stats]
    @test scalar_timeseries(cold_spell_max_length(cold_spell_cube; window=3)) == [stats[3] for stats in expected_cold_spell_stats]
    @test scalar_timeseries(cold_spell_total_length(cold_spell_cube; window=3)) == [stats[2] for stats in expected_cold_spell_stats]
    @test scalar_timeseries(hot_spell_frequency(hot_spell_cube; window=3)) == [stats[1] for stats in expected_hot_spell_stats]
    @test scalar_timeseries(hot_spell_max_length(hot_spell_cube; window=3)) == [stats[3] for stats in expected_hot_spell_stats]
    @test scalar_timeseries(hot_spell_total_length(hot_spell_cube; window=3)) == [stats[2] for stats in expected_hot_spell_stats]
    @test scalar_timeseries(heat_wave_index(hot_spell_cube)) == expected_heat_wave_index

    heatwave_tasmin_pattern_2000 = [21.0, 23.0, 24.0, 25.0, 18.0, 23.0, 24.0, 25.0, 19.0, 20.0]
    heatwave_tasmax_pattern_2000 = [29.0, 31.0, 32.0, 33.0, 20.0, 34.0, 35.0, 36.0, 21.0, 22.0]
    heatwave_tasmin_pattern_2001 = [23.0, 24.0, 25.0, 26.0, 27.0, 18.0, 22.0, 23.0, 24.0, 25.0]
    heatwave_tasmax_pattern_2001 = [31.0, 32.0, 33.0, 34.0, 35.0, 20.0, 30.0, 31.0, 32.0, 33.0]

    heatwave_tasmin_values = Float64[]
    heatwave_tasmax_values = Float64[]
    for date in dates
        idx = mod(dayofyear(date) - 1, 10) + 1
        if year(date) == 2000
            push!(heatwave_tasmin_values, heatwave_tasmin_pattern_2000[idx])
            push!(heatwave_tasmax_values, heatwave_tasmax_pattern_2000[idx])
        else
            push!(heatwave_tasmin_values, heatwave_tasmin_pattern_2001[idx])
            push!(heatwave_tasmax_values, heatwave_tasmax_pattern_2001[idx])
        end
    end

    heatwave_tasmin = make_cube(heatwave_tasmin_values, dates)
    heatwave_tasmax = make_cube(heatwave_tasmax_values, dates)

    expected_tx_tn_days_above = [count_pair_matches(first_group, second_group, (first_value, second_value) -> first_value > 22.0 && second_value > 30.0) for (first_group, second_group) in zip(grouped(heatwave_tasmin_values, dates, "YS"), grouped(heatwave_tasmax_values, dates, "YS"))]
    expected_heatwave_stats = [pair_run_length_stats(first_group, second_group, (first_value, second_value) -> first_value > 22.0 && second_value > 30.0, 3) for (first_group, second_group) in zip(grouped(heatwave_tasmin_values, dates, "YS"), grouped(heatwave_tasmax_values, dates, "YS"))]

    @test scalar_timeseries(tx_tn_days_above(heatwave_tasmin, heatwave_tasmax)) == expected_tx_tn_days_above
    @test scalar_timeseries(heat_wave_frequency(heatwave_tasmin, heatwave_tasmax)) == [stats[1] for stats in expected_heatwave_stats]
    @test scalar_timeseries(heat_wave_max_length(heatwave_tasmin, heatwave_tasmax)) == [stats[3] for stats in expected_heatwave_stats]
    @test scalar_timeseries(heat_wave_total_length(heatwave_tasmin, heatwave_tasmax)) == [stats[2] for stats in expected_heatwave_stats]

    high_precip_pattern_2000 = [0.0, 0.5, 1.0, 2.0, 0.3, 1.5, 2.5, 0.1, 1.0, 3.0]
    low_temp_pattern_2000 = [0.0, -1.0, -2.0, 1.0, -0.5, -1.5, 2.0, -3.0, -0.1, -4.0]
    high_precip_pattern_2001 = [0.2, 1.0, 2.0, 0.3, 1.2, 0.0, 3.0, 1.1, 0.9, 2.5]
    low_temp_pattern_2001 = [-1.0, -0.3, 0.5, -2.0, -1.0, -3.0, -0.1, 1.0, -0.5, -4.0]

    high_precip_values = Float64[]
    low_temp_values = Float64[]
    for date in dates
        idx = mod(dayofyear(date) - 1, 10) + 1
        if year(date) == 2000
            push!(high_precip_values, high_precip_pattern_2000[idx])
            push!(low_temp_values, low_temp_pattern_2000[idx])
        else
            push!(high_precip_values, high_precip_pattern_2001[idx])
            push!(low_temp_values, low_temp_pattern_2001[idx])
        end
    end

    high_precip_cube = make_cube(high_precip_values, dates)
    low_temp_cube = make_cube(low_temp_values, dates)
    expected_high_precip_low_temp = [count_pair_matches(pr_group, tas_group, (pr_value, tas_value) -> pr_value >= 1.0 && tas_value < 0.0) for (pr_group, tas_group) in zip(grouped(high_precip_values, dates, "YS"), grouped(low_temp_values, dates, "YS"))]

    @test scalar_timeseries(high_precip_low_temp(high_precip_cube, low_temp_cube; pr_thresh=1.0, tas_thresh=0.0)) == expected_high_precip_low_temp

    percentile_tx_pattern_2000 = [5.0, 9.0, 10.0, 11.0, 12.0, 6.0, 7.0, 13.0, 14.0, 5.0]
    percentile_tx_pattern_2001 = [4.0, 5.0, 6.0, 15.0, 16.0, 17.0, 3.0, 2.0, 18.0, 19.0]
    percentile_tn_pattern_2000 = [2.0, 1.0, 0.0, 5.0, 6.0, 2.0, 1.0, 7.0, 8.0, 9.0]
    percentile_tn_pattern_2001 = [3.0, 2.0, 1.0, 6.0, 7.0, 8.0, 2.0, 1.0, 9.0, 10.0]
    percentile_tg_pattern_2000 = [4.0, 5.0, 6.0, 7.0, 8.0, 4.0, 5.0, 8.0, 9.0, 3.0]
    percentile_tg_pattern_2001 = [2.0, 3.0, 4.0, 9.0, 10.0, 11.0, 3.0, 2.0, 12.0, 13.0]

    tx90_threshold_pattern_2000 = [8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 12.0, 12.0, 8.0]
    tx90_threshold_pattern_2001 = [7.0, 7.0, 7.0, 12.0, 12.0, 12.0, 7.0, 7.0, 17.0, 17.0]
    tx10_threshold_pattern_2000 = [6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 10.0, 10.0, 6.0]
    tx10_threshold_pattern_2001 = [4.5, 4.5, 4.5, 13.0, 13.0, 13.0, 4.5, 4.5, 16.0, 16.0]

    tn90_threshold_pattern_2000 = [1.0, 1.0, -1.0, 4.0, 4.0, 1.0, 1.0, 6.0, 6.0, 7.0]
    tn90_threshold_pattern_2001 = [2.0, 2.0, 0.0, 5.0, 5.0, 6.0, 1.0, 1.0, 8.0, 8.0]
    tn10_threshold_pattern_2000 = [3.0, 2.0, 1.0, 4.0, 5.0, 3.0, 2.0, 6.0, 7.0, 8.0]
    tn10_threshold_pattern_2001 = [4.0, 3.0, 2.0, 5.0, 6.0, 7.0, 3.0, 2.0, 8.0, 9.0]

    tg90_threshold_pattern_2000 = [3.0, 4.0, 5.0, 6.0, 6.0, 3.0, 4.0, 7.0, 8.0, 3.0]
    tg90_threshold_pattern_2001 = [1.0, 2.0, 3.0, 8.0, 9.0, 9.0, 2.0, 1.0, 10.0, 11.0]
    tg10_threshold_pattern_2000 = [4.5, 5.5, 6.5, 6.5, 7.5, 4.5, 5.5, 7.5, 8.5, 3.5]
    tg10_threshold_pattern_2001 = [2.5, 3.5, 4.5, 8.5, 9.5, 10.5, 3.5, 2.5, 11.5, 12.5]

    percentile_tx_values = Float64[]
    percentile_tn_values = Float64[]
    percentile_tg_values = Float64[]
    tx90_threshold_values = Float64[]
    tx10_threshold_values = Float64[]
    tn90_threshold_values = Float64[]
    tn10_threshold_values = Float64[]
    tg90_threshold_values = Float64[]
    tg10_threshold_values = Float64[]

    for date in dates
        idx = mod(dayofyear(date) - 1, 10) + 1
        if year(date) == 2000
            push!(percentile_tx_values, percentile_tx_pattern_2000[idx])
            push!(percentile_tn_values, percentile_tn_pattern_2000[idx])
            push!(percentile_tg_values, percentile_tg_pattern_2000[idx])
            push!(tx90_threshold_values, tx90_threshold_pattern_2000[idx])
            push!(tx10_threshold_values, tx10_threshold_pattern_2000[idx])
            push!(tn90_threshold_values, tn90_threshold_pattern_2000[idx])
            push!(tn10_threshold_values, tn10_threshold_pattern_2000[idx])
            push!(tg90_threshold_values, tg90_threshold_pattern_2000[idx])
            push!(tg10_threshold_values, tg10_threshold_pattern_2000[idx])
        else
            push!(percentile_tx_values, percentile_tx_pattern_2001[idx])
            push!(percentile_tn_values, percentile_tn_pattern_2001[idx])
            push!(percentile_tg_values, percentile_tg_pattern_2001[idx])
            push!(tx90_threshold_values, tx90_threshold_pattern_2001[idx])
            push!(tx10_threshold_values, tx10_threshold_pattern_2001[idx])
            push!(tn90_threshold_values, tn90_threshold_pattern_2001[idx])
            push!(tn10_threshold_values, tn10_threshold_pattern_2001[idx])
            push!(tg90_threshold_values, tg90_threshold_pattern_2001[idx])
            push!(tg10_threshold_values, tg10_threshold_pattern_2001[idx])
        end
    end

    percentile_tx = make_cube(percentile_tx_values, dates)
    percentile_tn = make_cube(percentile_tn_values, dates)
    percentile_tg = make_cube(percentile_tg_values, dates)
    tx90_threshold = make_cube(tx90_threshold_values, dates)
    tx10_threshold = make_cube(tx10_threshold_values, dates)
    tn90_threshold = make_cube(tn90_threshold_values, dates)
    tn10_threshold = make_cube(tn10_threshold_values, dates)
    tg90_threshold = make_cube(tg90_threshold_values, dates)
    tg10_threshold = make_cube(tg10_threshold_values, dates)

    expected_tx90p = [count_matches_against(group, thresh, >) for (group, thresh) in zip(grouped(percentile_tx_values, dates, "YS"), grouped(tx90_threshold_values, dates, "YS"))]
    expected_tx10p = [count_matches_against(group, thresh, <) for (group, thresh) in zip(grouped(percentile_tx_values, dates, "YS"), grouped(tx10_threshold_values, dates, "YS"))]
    expected_tn90p = [count_matches_against(group, thresh, >) for (group, thresh) in zip(grouped(percentile_tn_values, dates, "YS"), grouped(tn90_threshold_values, dates, "YS"))]
    expected_tn10p = [count_matches_against(group, thresh, <) for (group, thresh) in zip(grouped(percentile_tn_values, dates, "YS"), grouped(tn10_threshold_values, dates, "YS"))]
    expected_tg90p = [count_matches_against(group, thresh, >) for (group, thresh) in zip(grouped(percentile_tg_values, dates, "YS"), grouped(tg90_threshold_values, dates, "YS"))]
    expected_tg10p = [count_matches_against(group, thresh, <) for (group, thresh) in zip(grouped(percentile_tg_values, dates, "YS"), grouped(tg10_threshold_values, dates, "YS"))]
    expected_csdi = [spell_total_length_against(group, thresh, <, 3) for (group, thresh) in zip(grouped(percentile_tn_values, dates, "YS"), grouped(tn10_threshold_values, dates, "YS"))]
    expected_wsdi = [spell_total_length_against(group, thresh, >, 3) for (group, thresh) in zip(grouped(percentile_tx_values, dates, "YS"), grouped(tx90_threshold_values, dates, "YS"))]

    @test scalar_timeseries(tx90p(percentile_tx, tx90_threshold)) == expected_tx90p
    @test scalar_timeseries(tx10p(percentile_tx, tx10_threshold)) == expected_tx10p
    @test scalar_timeseries(tn90p(percentile_tn, tn90_threshold)) == expected_tn90p
    @test scalar_timeseries(tn10p(percentile_tn, tn10_threshold)) == expected_tn10p
    @test scalar_timeseries(tg90p(percentile_tg, tg90_threshold)) == expected_tg90p
    @test scalar_timeseries(tg10p(percentile_tg, tg10_threshold)) == expected_tg10p
    @test scalar_timeseries(cold_spell_duration_index(percentile_tn, tn10_threshold; window=3)) == expected_csdi
    @test scalar_timeseries(warm_spell_duration_index(percentile_tx, tx90_threshold; window=3)) == expected_wsdi
    @test_throws ErrorException tx90p(short_cube, percentile_tx)
end