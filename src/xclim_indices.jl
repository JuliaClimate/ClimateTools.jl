function _normalize_resample_frequency(freq)
    normalized = lowercase(string(freq))
    if normalized in ("ys", "as", "year", "yearly")
        return :yearly
    elseif normalized in ("ms", "month", "monthly")
        return :monthly
    end

    error("Unsupported freq=$(freq). Supported values are \"YS\" and \"MS\".")
end

function _period_groups(time_values, freq)
    frequency = _normalize_resample_frequency(freq)

    if frequency == :yearly
        keys = year.(time_values)
        unique_keys = unique(keys)
        dates = dates_builder_yearmonthday_hardcode(unique_keys, imois=7, iday=1)
    else
        keys = yearmonth.(time_values)
        unique_keys = unique(keys)
        dates = dates_builder_monthly_resample(unique_keys)
    end

    index_list = [findall(==(key), keys) for key in unique_keys]
    return index_list, dates
end

function _reduce_valid(data, reducer::Function)
    values = _valid_values(data)
    return isempty(values) ? NaN : reducer(values)
end

function _period_reduce_kernel(xout, xin; index_list, reducer::Function)
    xout .= NaN

    for i in eachindex(index_list)
        value = reducer(view(xin, index_list[i]))
        if !(ismissing(value) || (value isa Number && isnan(value)))
            xout[i] = value
        end
    end
end

function _period_reduce(cube::YAXArray; reducer::Function, freq="YS")
    index_list, dates = _period_groups(cube.time, freq)
    return _xmap_call(
        _period_reduce_kernel,
        cube;
        reduced_dims=:time,
        output_axes=(Dim{:time}(dates),),
        function_kwargs=(index_list=index_list, reducer=reducer),
    )
end

function _rolling_sum_kernel(xout, xin; window::Int)
    xout .= NaN

    for i in eachindex(xin)
        first_index = i - window + 1
        if first_index < 1
            continue
        end

        values = _valid_values(view(xin, first_index:i))
        if length(values) == window
            xout[i] = sum(values)
        end
    end
end

function _rolling_sum(cube::YAXArray; window::Int)
    return _xmap_call(
        _rolling_sum_kernel,
        cube;
        reduced_dims=:time,
        output_axes=(Dim{:time}(lookup(cube, :time)),),
        function_kwargs=(window=window,),
    )
end

function _resolve_reducer(op::Function)
    return op
end

function _resolve_reducer(op)
    normalized = lowercase(string(op))
    if normalized == "min"
        return Base.minimum
    elseif normalized == "max"
        return Base.maximum
    elseif normalized == "mean"
        return mean
    elseif normalized == "std"
        return std
    end

    error("Unsupported reducer $(op). Supported values are min, max, mean, std, or a function.")
end

function _comparison_operator(op)
    normalized = lowercase(string(op))
    if normalized in (">", "gt")
        return (value, threshold) -> value > threshold
    elseif normalized in (">=", "ge")
        return (value, threshold) -> value >= threshold
    elseif normalized in ("<", "lt")
        return (value, threshold) -> value < threshold
    elseif normalized in ("<=", "le")
        return (value, threshold) -> value <= threshold
    end

    error("Unsupported comparison operator $(op).")
end

function _threshold_count(cube::YAXArray; thresh, freq="YS", op=">=")
    comparator = _comparison_operator(op)
    reducer = data -> begin
        values = _valid_values(data)
        return isempty(values) ? NaN : count(value -> comparator(value, thresh), values)
    end

    return _period_reduce(cube; freq=freq, reducer=reducer)
end

function _threshold_proportion(cube::YAXArray; thresh, freq="YS", op=">=")
    comparator = _comparison_operator(op)
    reducer = data -> begin
        values = _valid_values(data)
        if isempty(values)
            return NaN
        end
        return count(value -> comparator(value, thresh), values) / length(values)
    end

    return _period_reduce(cube; freq=freq, reducer=reducer)
end

function _longest_run_length(data, comparator::Function, thresh)
    best = 0
    current = 0
    has_valid = false

    for value in data
        if ismissing(value) || (value isa Number && isnan(value))
            current = 0
            continue
        end

        has_valid = true
        if comparator(value, thresh)
            current += 1
            best = max(best, current)
        else
            current = 0
        end
    end

    return has_valid ? best : NaN
end

function _longest_threshold_run(cube::YAXArray; thresh, freq="YS", op)
    comparator = _comparison_operator(op)
    reducer = data -> _longest_run_length(data, comparator, thresh)
    return _period_reduce(cube; freq=freq, reducer=reducer)
end

function _thresholded_sum(cube::YAXArray; thresh, freq="YS", op=">")
    comparator = _comparison_operator(op)
    normalized = lowercase(string(op))
    difference = if normalized in (">", ">=", "gt", "ge")
        (value, threshold) -> value - threshold
    else
        (value, threshold) -> threshold - value
    end

    reducer = data -> begin
        total = 0.0
        has_valid = false

        for value in data
            if ismissing(value) || (value isa Number && isnan(value))
                continue
            end

            has_valid = true
            if comparator(value, thresh)
                total += difference(value, thresh)
            end
        end

        return has_valid ? total : NaN
    end

    return _period_reduce(cube; freq=freq, reducer=reducer)
end

function _spell_length_stats(data, comparator::Function, thresh, window::Int)
    window >= 1 || error("window must be >= 1")

    event_count = 0
    total_days = 0
    max_length = 0
    current = 0
    has_valid = false

    for value in data
        if ismissing(value) || (value isa Number && isnan(value))
            if current >= window
                event_count += 1
                total_days += current
                max_length = max(max_length, current)
            end
            current = 0
            continue
        end

        has_valid = true
        if comparator(value, thresh)
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

function _spell_stat_threshold(cube::YAXArray; thresh, freq="YS", op, window::Int=1, stat::Symbol)
    comparator = _comparison_operator(op)
    reducer = data -> begin
        event_count, total_days, max_length = _spell_length_stats(data, comparator, thresh, window)
        if stat == :count
            return event_count
        elseif stat == :sum
            return total_days
        elseif stat == :max
            return max_length
        end

        error("Unsupported spell statistic $(stat).")
    end

    return _period_reduce(cube; freq=freq, reducer=reducer)
end

function _threshold_count_pair_kernel(xout, xfirst, xsecond; index_list, comparator_first::Function, threshold_first, comparator_second::Function, threshold_second)
    xout .= NaN

    for i in eachindex(index_list)
        matches = 0
        valid = 0

        for idx in index_list[i]
            first_value = xfirst[idx]
            second_value = xsecond[idx]

            if ismissing(first_value) || ismissing(second_value) ||
               (first_value isa Number && isnan(first_value)) ||
               (second_value isa Number && isnan(second_value))
                continue
            end

            valid += 1
            if comparator_first(first_value, threshold_first) && comparator_second(second_value, threshold_second)
                matches += 1
            end
        end

        if valid > 0
            xout[i] = matches
        end
    end
end

function _threshold_count_pair(first::YAXArray, second::YAXArray; first_thresh, second_thresh, freq="YS", first_op, second_op)
    _ensure_matching_time_axis(first, second)
    index_list, dates = _period_groups(first.time, freq)
    return _xmap_call(
        _threshold_count_pair_kernel,
        (first, second);
        reduced_dims=:time,
        output_axes=(Dim{:time}(dates),),
        function_kwargs=(
            index_list=index_list,
            comparator_first=_comparison_operator(first_op),
            threshold_first=first_thresh,
            comparator_second=_comparison_operator(second_op),
            threshold_second=second_thresh,
        ),
    )
end

function _spell_length_stats_pair(first_data, second_data, comparator_first::Function, thresh_first, comparator_second::Function, thresh_second, window::Int)
    window >= 1 || error("window must be >= 1")

    event_count = 0
    total_days = 0
    max_length = 0
    current = 0
    has_valid = false

    for (first_value, second_value) in zip(first_data, second_data)
        if ismissing(first_value) || ismissing(second_value) ||
           (first_value isa Number && isnan(first_value)) ||
           (second_value isa Number && isnan(second_value))
            if current >= window
                event_count += 1
                total_days += current
                max_length = max(max_length, current)
            end
            current = 0
            continue
        end

        has_valid = true
        if comparator_first(first_value, thresh_first) && comparator_second(second_value, thresh_second)
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

function _spell_stat_pair_kernel(xout, xfirst, xsecond; index_list, comparator_first::Function, threshold_first, comparator_second::Function, threshold_second, window::Int, stat::Symbol)
    xout .= NaN

    for i in eachindex(index_list)
        event_count, total_days, max_length = _spell_length_stats_pair(
            view(xfirst, index_list[i]),
            view(xsecond, index_list[i]),
            comparator_first,
            threshold_first,
            comparator_second,
            threshold_second,
            window,
        )

        if stat == :count
            xout[i] = event_count
        elseif stat == :sum
            xout[i] = total_days
        elseif stat == :max
            xout[i] = max_length
        else
            error("Unsupported spell statistic $(stat).")
        end
    end
end

function _spell_stat_pair(first::YAXArray, second::YAXArray; first_thresh, second_thresh, freq="YS", first_op, second_op, window::Int=1, stat::Symbol)
    _ensure_matching_time_axis(first, second)
    index_list, dates = _period_groups(first.time, freq)
    return _xmap_call(
        _spell_stat_pair_kernel,
        (first, second);
        reduced_dims=:time,
        output_axes=(Dim{:time}(dates),),
        function_kwargs=(
            index_list=index_list,
            comparator_first=_comparison_operator(first_op),
            threshold_first=first_thresh,
            comparator_second=_comparison_operator(second_op),
            threshold_second=second_thresh,
            window=window,
            stat=stat,
        ),
    )
end

function _threshold_count_against_kernel(xout, xin, xthreshold; index_list, comparator::Function)
    xout .= NaN

    for i in eachindex(index_list)
        matches = 0
        valid = 0

        for idx in index_list[i]
            value = xin[idx]
            threshold = xthreshold[idx]

            if ismissing(value) || ismissing(threshold) ||
               (value isa Number && isnan(value)) ||
               (threshold isa Number && isnan(threshold))
                continue
            end

            valid += 1
            if comparator(value, threshold)
                matches += 1
            end
        end

        if valid > 0
            xout[i] = matches
        end
    end
end

function _windowed_run_count_against(data, thresholds, comparator::Function, window::Int)
    total = 0
    current = 0
    has_valid = false

    for (value, threshold) in zip(data, thresholds)
        if ismissing(value) || ismissing(threshold) ||
           (value isa Number && isnan(value)) ||
           (threshold isa Number && isnan(threshold))
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

function _windowed_threshold_run_count_against_kernel(xout, xin, xthreshold; index_list, comparator::Function, window::Int)
    xout .= NaN

    for i in eachindex(index_list)
        xout[i] = _windowed_run_count_against(view(xin, index_list[i]), view(xthreshold, index_list[i]), comparator, window)
    end
end

function _ensure_matching_size(first::YAXArray, second::YAXArray)
    size(first) == size(second) || error("Input cubes must have the same size.")
end

function _ensure_matching_time_axis(first::YAXArray, second::YAXArray)
    _ensure_matching_size(first, second)
    collect(first.time) == collect(second.time) || error("Input cubes must share the same time axis.")
end

function _threshold_count_against(cube::YAXArray, threshold_cube::YAXArray; freq="YS", op=">")
    _ensure_matching_time_axis(cube, threshold_cube)
    index_list, dates = _period_groups(cube.time, freq)
    return _xmap_call(
        _threshold_count_against_kernel,
        (cube, threshold_cube);
        reduced_dims=:time,
        output_axes=(Dim{:time}(dates),),
        function_kwargs=(index_list=index_list, comparator=_comparison_operator(op)),
    )
end

function _windowed_threshold_run_count_against(cube::YAXArray, threshold_cube::YAXArray; window::Int=6, freq="YS", op=">")
    window >= 1 || error("window must be >= 1")
    _ensure_matching_time_axis(cube, threshold_cube)
    index_list, dates = _period_groups(cube.time, freq)
    return _xmap_call(
        _windowed_threshold_run_count_against_kernel,
        (cube, threshold_cube);
        reduced_dims=:time,
        output_axes=(Dim{:time}(dates),),
        function_kwargs=(index_list=index_list, comparator=_comparison_operator(op), window=window),
    )
end

function _temperature_range_cube(tasmin::YAXArray, tasmax::YAXArray)
    _ensure_matching_size(tasmin, tasmax)
    return _cube_with_data(tasmax, Array(tasmax) .- Array(tasmin))
end

"""
    tg_max(tas::YAXArray; freq="YS")

Highest daily mean temperature over each resampling period.
"""
function tg_max(tas::YAXArray; freq="YS")
    return _period_reduce(tas; freq=freq, reducer=data -> _reduce_valid(data, Base.maximum))
end

"""
    tg_mean(tas::YAXArray; freq="YS")

Mean daily mean temperature over each resampling period.
"""
function tg_mean(tas::YAXArray; freq="YS")
    return _period_reduce(tas; freq=freq, reducer=data -> _reduce_valid(data, mean))
end

"""
    tg_min(tas::YAXArray; freq="YS")

Lowest daily mean temperature over each resampling period.
"""
function tg_min(tas::YAXArray; freq="YS")
    return _period_reduce(tas; freq=freq, reducer=data -> _reduce_valid(data, Base.minimum))
end

"""
    tx_max(tasmax::YAXArray; freq="YS")

Highest daily maximum temperature over each resampling period.
"""
function tx_max(tasmax::YAXArray; freq="YS")
    return _period_reduce(tasmax; freq=freq, reducer=data -> _reduce_valid(data, Base.maximum))
end

"""
    tx_min(tasmax::YAXArray; freq="YS")

Lowest daily maximum temperature over each resampling period.
"""
function tx_min(tasmax::YAXArray; freq="YS")
    return _period_reduce(tasmax; freq=freq, reducer=data -> _reduce_valid(data, Base.minimum))
end

"""
    tn_max(tasmin::YAXArray; freq="YS")

Highest daily minimum temperature over each resampling period.
"""
function tn_max(tasmin::YAXArray; freq="YS")
    return _period_reduce(tasmin; freq=freq, reducer=data -> _reduce_valid(data, Base.maximum))
end

"""
    tn_min(tasmin::YAXArray; freq="YS")

Lowest daily minimum temperature over each resampling period.
"""
function tn_min(tasmin::YAXArray; freq="YS")
    return _period_reduce(tasmin; freq=freq, reducer=data -> _reduce_valid(data, Base.minimum))
end

"""
    daily_temperature_range(tasmin::YAXArray, tasmax::YAXArray; freq="YS", op="mean")

Statistic of the daily temperature range over each resampling period.
"""
function daily_temperature_range(tasmin::YAXArray, tasmax::YAXArray; freq="YS", op="mean")
    reducer = _resolve_reducer(op)
    dtr = _temperature_range_cube(tasmin, tasmax)
    return _period_reduce(dtr; freq=freq, reducer=data -> _reduce_valid(data, reducer))
end

"""
    daily_temperature_range_variability(tasmin::YAXArray, tasmax::YAXArray; freq="YS")

Mean absolute day-to-day variation in daily temperature range.
"""
function daily_temperature_range_variability(tasmin::YAXArray, tasmax::YAXArray; freq="YS")
    dtr = _temperature_range_cube(tasmin, tasmax)
    reducer = data -> begin
        values = _valid_values(data)
        if length(values) < 2
            return NaN
        end
        return mean(abs.(Base.diff(values)))
    end

    return _period_reduce(dtr; freq=freq, reducer=reducer)
end

"""
    extreme_temperature_range(tasmin::YAXArray, tasmax::YAXArray; freq="YS")

Maximum daily maximum temperature minus minimum daily minimum temperature for each period.
"""
function extreme_temperature_range(tasmin::YAXArray, tasmax::YAXArray; freq="YS")
    txx = tx_max(tasmax; freq=freq)
    tnn = tn_min(tasmin; freq=freq)
    return _cube_with_data(txx, Array(txx) .- Array(tnn))
end

"""
    max_1day_precipitation_amount(pr::YAXArray; freq="YS")

Highest 1-day precipitation amount over each resampling period.
"""
function max_1day_precipitation_amount(pr::YAXArray; freq="YS")
    return _period_reduce(pr; freq=freq, reducer=data -> _reduce_valid(data, Base.maximum))
end

"""
    max_n_day_precipitation_amount(pr::YAXArray; window=1, freq="YS")

Highest rolling n-day precipitation amount over each resampling period.
"""
function max_n_day_precipitation_amount(pr::YAXArray; window::Int=1, freq="YS")
    window >= 1 || error("window must be >= 1")
    rolling = _rolling_sum(pr; window=window)
    return _period_reduce(rolling; freq=freq, reducer=data -> _reduce_valid(data, Base.maximum))
end

"""
    daily_pr_intensity(pr::YAXArray; thresh=1.0, freq="YS", op=">=")

Average precipitation over wet days for each resampling period.
"""
function daily_pr_intensity(pr::YAXArray; thresh=1.0, freq="YS", op=">=")
    comparator = _comparison_operator(op)
    reducer = data -> begin
        values = _valid_values(data)
        wet_values = [value for value in values if comparator(value, thresh)]
        return isempty(wet_values) ? NaN : mean(wet_values)
    end

    return _period_reduce(pr; freq=freq, reducer=reducer)
end

"""
    hot_days(tasmax::YAXArray; thresh=25.0, freq="YS", op=">")

Number of days where daily maximum temperature exceeds a threshold.
"""
function hot_days(tasmax::YAXArray; thresh=25.0, freq="YS", op=">")
    return _threshold_count(tasmax; thresh=thresh, freq=freq, op=op)
end

"""
    frost_days(tasmin::YAXArray; thresh=0.0, freq="YS", op="<")

Number of days where daily minimum temperature is below a frost threshold.
"""
function frost_days(tasmin::YAXArray; thresh=0.0, freq="YS", op="<")
    return _threshold_count(tasmin; thresh=thresh, freq=freq, op=op)
end

"""
    tg_days_above(tas::YAXArray; thresh=10.0, freq="YS", op=">")

Number of days where daily mean temperature exceeds a threshold.
"""
function tg_days_above(tas::YAXArray; thresh=10.0, freq="YS", op=">")
    return _threshold_count(tas; thresh=thresh, freq=freq, op=op)
end

"""
    tg_days_below(tas::YAXArray; thresh=10.0, freq="YS", op="<")

Number of days where daily mean temperature is below a threshold.
"""
function tg_days_below(tas::YAXArray; thresh=10.0, freq="YS", op="<")
    return _threshold_count(tas; thresh=thresh, freq=freq, op=op)
end

"""
    tn_days_above(tasmin::YAXArray; thresh=20.0, freq="YS", op=">")

Number of days where daily minimum temperature exceeds a threshold.
"""
function tn_days_above(tasmin::YAXArray; thresh=20.0, freq="YS", op=">")
    return _threshold_count(tasmin; thresh=thresh, freq=freq, op=op)
end

"""
    tn_days_below(tasmin::YAXArray; thresh=-10.0, freq="YS", op="<")

Number of days where daily minimum temperature is below a threshold.
"""
function tn_days_below(tasmin::YAXArray; thresh=-10.0, freq="YS", op="<")
    return _threshold_count(tasmin; thresh=thresh, freq=freq, op=op)
end

"""
    tx_days_above(tasmax::YAXArray; thresh=25.0, freq="YS", op=">")

Number of days where daily maximum temperature exceeds a threshold.
"""
function tx_days_above(tasmax::YAXArray; thresh=25.0, freq="YS", op=">")
    return _threshold_count(tasmax; thresh=thresh, freq=freq, op=op)
end

"""
    tx_days_below(tasmax::YAXArray; thresh=25.0, freq="YS", op="<")

Number of days where daily maximum temperature is below a threshold.
"""
function tx_days_below(tasmax::YAXArray; thresh=25.0, freq="YS", op="<")
    return _threshold_count(tasmax; thresh=thresh, freq=freq, op=op)
end

"""
    warm_day_frequency(tasmax::YAXArray; thresh=30.0, freq="YS", op=">")

Number of days where daily maximum temperature exceeds a warm-day threshold.
"""
function warm_day_frequency(tasmax::YAXArray; thresh=30.0, freq="YS", op=">")
    return _threshold_count(tasmax; thresh=thresh, freq=freq, op=op)
end

"""
    warm_night_frequency(tasmin::YAXArray; thresh=22.0, freq="YS", op=">")

Number of days where daily minimum temperature exceeds a warm-night threshold.
"""
function warm_night_frequency(tasmin::YAXArray; thresh=22.0, freq="YS", op=">")
    return _threshold_count(tasmin; thresh=thresh, freq=freq, op=op)
end

"""
    ice_days(tasmax::YAXArray; thresh=0.0, freq="YS", op="<")

Number of days where daily maximum temperature remains below freezing.
"""
function ice_days(tasmax::YAXArray; thresh=0.0, freq="YS", op="<")
    return _threshold_count(tasmax; thresh=thresh, freq=freq, op=op)
end

"""
    growing_degree_days(tas::YAXArray; thresh=4.0, freq="YS")

Sum of positive daily mean temperature exceedances above a growing-degree threshold.
"""
function growing_degree_days(tas::YAXArray; thresh=4.0, freq="YS")
    return _thresholded_sum(tas; thresh=thresh, freq=freq, op=">")
end

"""
    heating_degree_days(tas::YAXArray; thresh=17.0, freq="YS")

Sum of temperature deficits below a heating threshold.
"""
function heating_degree_days(tas::YAXArray; thresh=17.0, freq="YS")
    return _thresholded_sum(tas; thresh=thresh, freq=freq, op="<")
end

"""
    cooling_degree_days(tas::YAXArray; thresh=18.0, freq="YS")

Sum of positive daily mean temperature exceedances above a cooling threshold.
"""
function cooling_degree_days(tas::YAXArray; thresh=18.0, freq="YS")
    return _thresholded_sum(tas; thresh=thresh, freq=freq, op=">")
end

"""
    dry_days(pr::YAXArray; thresh=0.2, freq="YS", op="<")

Number of days where daily precipitation is below a threshold.
"""
function dry_days(pr::YAXArray; thresh=0.2, freq="YS", op="<")
    return _threshold_count(pr; thresh=thresh, freq=freq, op=op)
end

"""
    wetdays(pr::YAXArray; thresh=1.0, freq="YS", op=">=")

Number of days where daily precipitation meets or exceeds a wet-day threshold.
"""
function wetdays(pr::YAXArray; thresh=1.0, freq="YS", op=">=")
    return _threshold_count(pr; thresh=thresh, freq=freq, op=op)
end

"""
    wetdays_prop(pr::YAXArray; thresh=1.0, freq="YS", op=">=")

Proportion of valid days meeting a wet-day threshold within each period.
"""
function wetdays_prop(pr::YAXArray; thresh=1.0, freq="YS", op=">=")
    return _threshold_proportion(pr; thresh=thresh, freq=freq, op=op)
end

"""
    maximum_consecutive_dry_days(pr::YAXArray; thresh=1.0, freq="YS")

Longest run of days with precipitation below a dry-day threshold.
"""
function maximum_consecutive_dry_days(pr::YAXArray; thresh=1.0, freq="YS")
    return _longest_threshold_run(pr; thresh=thresh, freq=freq, op="<")
end

"""
    maximum_consecutive_wet_days(pr::YAXArray; thresh=1.0, freq="YS")

Longest run of days with precipitation above a wet-day threshold.
"""
function maximum_consecutive_wet_days(pr::YAXArray; thresh=1.0, freq="YS")
    return _longest_threshold_run(pr; thresh=thresh, freq=freq, op=">")
end

"""
    maximum_consecutive_frost_days(tasmin::YAXArray; thresh=0.0, freq="YS")

Longest run of days with daily minimum temperature below a frost threshold.
"""
function maximum_consecutive_frost_days(tasmin::YAXArray; thresh=0.0, freq="YS")
    return _longest_threshold_run(tasmin; thresh=thresh, freq=freq, op="<")
end

"""
    maximum_consecutive_frost_free_days(tasmin::YAXArray; thresh=0.0, freq="YS")

Longest run of days with daily minimum temperature at or above a threshold.
"""
function maximum_consecutive_frost_free_days(tasmin::YAXArray; thresh=0.0, freq="YS")
    return _longest_threshold_run(tasmin; thresh=thresh, freq=freq, op=">=")
end

"""
    maximum_consecutive_tx_days(tasmax::YAXArray; thresh=25.0, freq="YS")

Longest run of days with daily maximum temperature above a threshold.
"""
function maximum_consecutive_tx_days(tasmax::YAXArray; thresh=25.0, freq="YS")
    return _longest_threshold_run(tasmax; thresh=thresh, freq=freq, op=">")
end

"""
    cold_spell_days(tas::YAXArray; thresh=-10.0, window=5, freq="YS", op="<")

Total number of days belonging to cold spells defined by a fixed temperature threshold.
"""
function cold_spell_days(tas::YAXArray; thresh=-10.0, window::Int=5, freq="YS", op="<")
    return _spell_stat_threshold(tas; thresh=thresh, freq=freq, op=op, window=window, stat=:sum)
end

"""
    cold_spell_frequency(tas::YAXArray; thresh=-10.0, window=5, freq="YS", op="<")

Number of cold spell events defined by a fixed temperature threshold.
"""
function cold_spell_frequency(tas::YAXArray; thresh=-10.0, window::Int=5, freq="YS", op="<")
    return _spell_stat_threshold(tas; thresh=thresh, freq=freq, op=op, window=window, stat=:count)
end

"""
    cold_spell_max_length(tas::YAXArray; thresh=-10.0, window=1, freq="YS", op="<")

Longest cold spell length defined by a fixed temperature threshold.
"""
function cold_spell_max_length(tas::YAXArray; thresh=-10.0, window::Int=1, freq="YS", op="<")
    return _spell_stat_threshold(tas; thresh=thresh, freq=freq, op=op, window=window, stat=:max)
end

"""
    cold_spell_total_length(tas::YAXArray; thresh=-10.0, window=3, freq="YS", op="<")

Total number of days belonging to cold spells defined by a fixed temperature threshold.
"""
function cold_spell_total_length(tas::YAXArray; thresh=-10.0, window::Int=3, freq="YS", op="<")
    return _spell_stat_threshold(tas; thresh=thresh, freq=freq, op=op, window=window, stat=:sum)
end

"""
    hot_spell_frequency(tasmax::YAXArray; thresh=30.0, window=3, freq="YS", op=">")

Number of hot spell events defined by a fixed maximum-temperature threshold.
"""
function hot_spell_frequency(tasmax::YAXArray; thresh=30.0, window::Int=3, freq="YS", op=">")
    return _spell_stat_threshold(tasmax; thresh=thresh, freq=freq, op=op, window=window, stat=:count)
end

"""
    hot_spell_max_length(tasmax::YAXArray; thresh=30.0, window=1, freq="YS", op=">")

Longest hot spell length defined by a fixed maximum-temperature threshold.
"""
function hot_spell_max_length(tasmax::YAXArray; thresh=30.0, window::Int=1, freq="YS", op=">")
    return _spell_stat_threshold(tasmax; thresh=thresh, freq=freq, op=op, window=window, stat=:max)
end

"""
    hot_spell_total_length(tasmax::YAXArray; thresh=30.0, window=3, freq="YS", op=">")

Total number of days belonging to hot spells defined by a fixed maximum-temperature threshold.
"""
function hot_spell_total_length(tasmax::YAXArray; thresh=30.0, window::Int=3, freq="YS", op=">")
    return _spell_stat_threshold(tasmax; thresh=thresh, freq=freq, op=op, window=window, stat=:sum)
end

"""
    heat_wave_index(tasmax::YAXArray; thresh=25.0, window=5, freq="YS", op=">")

Total number of days belonging to heat waves defined on daily maximum temperature alone.
"""
function heat_wave_index(tasmax::YAXArray; thresh=25.0, window::Int=5, freq="YS", op=">")
    return hot_spell_total_length(tasmax; thresh=thresh, window=window, freq=freq, op=op)
end

"""
    tx_tn_days_above(tasmin::YAXArray, tasmax::YAXArray; thresh_tasmin=22.0, thresh_tasmax=30.0, freq="YS", op=">")

Number of days where daily minimum and maximum temperatures both exceed their thresholds.
"""
function tx_tn_days_above(tasmin::YAXArray, tasmax::YAXArray; thresh_tasmin=22.0, thresh_tasmax=30.0, freq="YS", op=">")
    return _threshold_count_pair(tasmin, tasmax; first_thresh=thresh_tasmin, second_thresh=thresh_tasmax, freq=freq, first_op=op, second_op=op)
end

"""
    heat_wave_frequency(tasmin::YAXArray, tasmax::YAXArray; thresh_tasmin=22.0, thresh_tasmax=30.0, window=3, freq="YS", op=">")

Number of heat-wave events where daily minimum and maximum temperatures both exceed thresholds.
"""
function heat_wave_frequency(tasmin::YAXArray, tasmax::YAXArray; thresh_tasmin=22.0, thresh_tasmax=30.0, window::Int=3, freq="YS", op=">")
    return _spell_stat_pair(tasmin, tasmax; first_thresh=thresh_tasmin, second_thresh=thresh_tasmax, freq=freq, first_op=op, second_op=op, window=window, stat=:count)
end

"""
    heat_wave_max_length(tasmin::YAXArray, tasmax::YAXArray; thresh_tasmin=22.0, thresh_tasmax=30.0, window=3, freq="YS", op=">")

Maximum heat-wave length where daily minimum and maximum temperatures both exceed thresholds.
"""
function heat_wave_max_length(tasmin::YAXArray, tasmax::YAXArray; thresh_tasmin=22.0, thresh_tasmax=30.0, window::Int=3, freq="YS", op=">")
    return _spell_stat_pair(tasmin, tasmax; first_thresh=thresh_tasmin, second_thresh=thresh_tasmax, freq=freq, first_op=op, second_op=op, window=window, stat=:max)
end

"""
    heat_wave_total_length(tasmin::YAXArray, tasmax::YAXArray; thresh_tasmin=22.0, thresh_tasmax=30.0, window=3, freq="YS", op=">")

Total number of days belonging to heat-wave events where daily minimum and maximum temperatures both exceed thresholds.
"""
function heat_wave_total_length(tasmin::YAXArray, tasmax::YAXArray; thresh_tasmin=22.0, thresh_tasmax=30.0, window::Int=3, freq="YS", op=">")
    return _spell_stat_pair(tasmin, tasmax; first_thresh=thresh_tasmin, second_thresh=thresh_tasmax, freq=freq, first_op=op, second_op=op, window=window, stat=:sum)
end

"""
    high_precip_low_temp(pr::YAXArray, tas::YAXArray; pr_thresh=0.4, tas_thresh=-0.2, freq="YS")

Number of days where precipitation meets a threshold while temperature stays below a threshold.
"""
function high_precip_low_temp(pr::YAXArray, tas::YAXArray; pr_thresh=0.4, tas_thresh=-0.2, freq="YS")
    return _threshold_count_pair(pr, tas; first_thresh=pr_thresh, second_thresh=tas_thresh, freq=freq, first_op=">=", second_op="<")
end

"""
    tg90p(tas::YAXArray, tas_per::YAXArray; freq="YS", op=">")

Count of days where daily mean temperature exceeds an aligned percentile threshold series.
"""
function tg90p(tas::YAXArray, tas_per::YAXArray; freq="YS", op=">")
    return _threshold_count_against(tas, tas_per; freq=freq, op=op)
end

"""
    tg10p(tas::YAXArray, tas_per::YAXArray; freq="YS", op="<")

Count of days where daily mean temperature is below an aligned percentile threshold series.
"""
function tg10p(tas::YAXArray, tas_per::YAXArray; freq="YS", op="<")
    return _threshold_count_against(tas, tas_per; freq=freq, op=op)
end

"""
    tn90p(tasmin::YAXArray, tasmin_per::YAXArray; freq="YS", op=">")

Count of days where daily minimum temperature exceeds an aligned percentile threshold series.
"""
function tn90p(tasmin::YAXArray, tasmin_per::YAXArray; freq="YS", op=">")
    return _threshold_count_against(tasmin, tasmin_per; freq=freq, op=op)
end

"""
    tn10p(tasmin::YAXArray, tasmin_per::YAXArray; freq="YS", op="<")

Count of days where daily minimum temperature is below an aligned percentile threshold series.
"""
function tn10p(tasmin::YAXArray, tasmin_per::YAXArray; freq="YS", op="<")
    return _threshold_count_against(tasmin, tasmin_per; freq=freq, op=op)
end

"""
    tx90p(tasmax::YAXArray, tasmax_per::YAXArray; freq="YS", op=">")

Count of days where daily maximum temperature exceeds an aligned percentile threshold series.
"""
function tx90p(tasmax::YAXArray, tasmax_per::YAXArray; freq="YS", op=">")
    return _threshold_count_against(tasmax, tasmax_per; freq=freq, op=op)
end

"""
    tx10p(tasmax::YAXArray, tasmax_per::YAXArray; freq="YS", op="<")

Count of days where daily maximum temperature is below an aligned percentile threshold series.
"""
function tx10p(tasmax::YAXArray, tasmax_per::YAXArray; freq="YS", op="<")
    return _threshold_count_against(tasmax, tasmax_per; freq=freq, op=op)
end

"""
    cold_spell_duration_index(tasmin::YAXArray, tasmin_per::YAXArray; window=6, freq="YS", op="<")

Total number of days belonging to cold spells, where a spell is a run of at least `window`
consecutive days below an aligned percentile threshold series.
"""
function cold_spell_duration_index(tasmin::YAXArray, tasmin_per::YAXArray; window::Int=6, freq="YS", op="<")
    return _windowed_threshold_run_count_against(tasmin, tasmin_per; window=window, freq=freq, op=op)
end

"""
    warm_spell_duration_index(tasmax::YAXArray, tasmax_per::YAXArray; window=6, freq="YS", op=">")

Total number of days belonging to warm spells, where a spell is a run of at least `window`
consecutive days above an aligned percentile threshold series.
"""
function warm_spell_duration_index(tasmax::YAXArray, tasmax_per::YAXArray; window::Int=6, freq="YS", op=">")
    return _windowed_threshold_run_count_against(tasmax, tasmax_per; window=window, freq=freq, op=op)
end