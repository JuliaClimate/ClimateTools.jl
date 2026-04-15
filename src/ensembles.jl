# Portions of this file are adapted from xclim (https://github.com/Ouranosinc/xclim),
# Copyright 2018-2023 Ouranos Inc. and contributors, and are distributed under
# the Apache License, Version 2.0. This ClimateTools.jl version rewrites and
# modifies that material for YAXArrays-based Julia workflows. See
# LICENSES/xclim-APACHE-2.0.txt and LICENSE.md for details.

const _ENSEMBLE_QUANTILE_PARAMS = Dict(
    "interpolated_inverted_cdf" => (0.0, 1.0),
    "hazen" => (0.5, 0.5),
    "weibull" => (0.0, 0.0),
    "linear" => (1.0, 1.0),
    "median_unbiased" => (1 / 3, 1 / 3),
    "normal_unbiased" => (3 / 8, 3 / 8),
)

function _axis_name(axis)
    return Symbol(name(axis))
end

function _axis_names(cube::YAXArray)
    return [_axis_name(axis) for axis in cube.axes]
end

function _has_axis(cube::YAXArray, dim)
    target = _normalize_xmap_dim(dim)
    return target in _axis_names(cube)
end

function _axis_position(cube::YAXArray, dim)
    target = _normalize_xmap_dim(dim)
    return findfirst(==(target), _axis_names(cube))
end

function _require_axis_position(cube::YAXArray, dim)
    position = _axis_position(cube, dim)
    isnothing(position) && error("Cube must have a $(dim) axis.")
    return position
end

function _required_member_count(min_members, total_members::Int)
    if isnothing(min_members)
        return total_members
    end

    required = Int(min_members)
    required >= 1 || error("min_members must be at least 1.")
    return required
end

function _resolve_realization_weights(weights, total_members::Int, realization_dim)
    isnothing(weights) && return nothing

    resolved = if weights isa YAXArray
        ndims(weights) == 1 || error("weights must be one-dimensional.")
        _require_axis_position(weights, realization_dim)
        Float64.(vec(Array(weights)))
    elseif weights isa AbstractVector
        Float64.(collect(weights))
    else
        error("weights must be an AbstractVector or a one-dimensional YAXArray.")
    end

    length(resolved) == total_members || error("weights length must match the realization axis length.")
    any(ismissing.(resolved)) && error("weights cannot contain missing values.")
    any(isnan.(resolved)) && error("weights cannot contain NaN values.")
    any(resolved .< 0) && error("weights must be non-negative.")
    return resolved
end

function _slice_axis(cube::YAXArray, axis_position::Int, index::Int)
    selectors = ntuple(i -> i == axis_position ? index : Colon(), ndims(cube))
    return cube[selectors...]
end

function _dataset_from_pairs(variable_pairs)
    return Dataset(; variable_pairs...)
end

function _reshape_output(values::AbstractVector, output_shape::Tuple)
    if isempty(output_shape)
        return reshape(copy(values), ())
    end
    return reshape(copy(values), output_shape)
end

function _value_is_valid(value)
    return !(ismissing(value) || isnan(value))
end

function _float_or_nan(value)
    return ismissing(value) ? NaN : Float64(value)
end

function _valid_numeric_positions(values)
    return [index for index in eachindex(values) if _value_is_valid(values[index])]
end

function _valid_numeric_values(values)
    positions = _valid_numeric_positions(values)
    return Float64[values[index] for index in positions], positions
end

function _weighted_mean_std(values::AbstractVector{<:Real}, weights::AbstractVector{<:Real})
    weight_sum = sum(weights)
    if iszero(weight_sum)
        return NaN, NaN
    end

    mean_value = sum(values .* weights) / weight_sum
    variance = sum(((values .- mean_value) .^ 2) .* weights) / weight_sum
    return mean_value, sqrt(max(variance, 0.0))
end

function _normalize_quantile_method(method)
    normalized = lowercase(String(method))
    haskey(_ENSEMBLE_QUANTILE_PARAMS, normalized) || error("Unsupported percentile method $(method).")
    return normalized
end

function _percentile_suffix(value)
    rounded = round(Float64(value); digits=6)
    if isapprox(rounded, round(Int, rounded); atol=1e-6)
        return string(round(Int, rounded))
    end
    return replace(string(rounded), "." => "_")
end

function _weighted_h(n::Float64, q::AbstractVector{<:Real}, method::AbstractString)
    quantiles = Float64.(q)
    h = if method == "linear"
        (n - 1) .* quantiles .+ 1
    elseif method == "interpolated_inverted_cdf"
        n .* quantiles
    elseif method == "hazen"
        n .* quantiles .+ 0.5
    elseif method == "weibull"
        (n + 1) .* quantiles
    elseif method == "median_unbiased"
        (n + 1 / 3) .* quantiles .+ 1 / 3
    elseif method == "normal_unbiased"
        (n + 1 / 4) .* quantiles .+ 3 / 8
    else
        error("Unsupported weighted percentile method $(method).")
    end

    return clamp.(h, 1, n)
end

function _weighted_quantiles_1d(values::AbstractVector{<:Real}, weights::AbstractVector{<:Real}, q::AbstractVector{<:Real}; method="linear")
    length(values) == length(weights) || error("weights and values must have the same length.")

    valid = .!isnan.(values)
    filtered_values = Float64.(values[valid])
    filtered_weights = Float64.(weights[valid])

    nonzero = filtered_weights .!= 0
    filtered_values = filtered_values[nonzero]
    filtered_weights = filtered_weights[nonzero]

    if isempty(filtered_values)
        return fill(NaN, length(q))
    end

    weight_sum = sum(filtered_weights)
    iszero(weight_sum) && return fill(NaN, length(q))

    effective_n = weight_sum^2 / sum(filtered_weights .^ 2)
    sorter = sortperm(filtered_values)
    filtered_values = filtered_values[sorter]
    filtered_weights = filtered_weights[sorter] ./ weight_sum

    weights_cumulative = vcat(0.0, cumsum(filtered_weights))
    interpolation_points = _weighted_h(effective_n, q, method)
    output = similar(Float64.(q))

    for (index, h) in pairs(interpolation_points)
        u = max.((h - 1) / effective_n, min.(h / effective_n, weights_cumulative))
        v = u .* effective_n .- h .+ 1
        output[index] = sum(filtered_values .* diff(v))
    end

    return output
end

function _ensemble_stats_kernel(xout, xin)
    xout .= NaN

    values, _ = _valid_numeric_values(xin)
    isempty(values) && return

    xout[1] = Statistics.maximum(values)
    xout[2] = Statistics.mean(values)
    xout[3] = Statistics.minimum(values)
end

"""
    ensemble_stats(xout, xin)

Compute the maximum, mean, and minimum values from the input array `xin` and store them in `xout`.
"""
function ensemble_stats(xout, xin)
    _ensemble_stats_kernel(xout, xin)
end

"""
    ensemble_stats(cube::YAXArray; dim="time")

Compute the maximum, mean, and minimum values from a YAXArray along dimension `dim`.
"""
function ensemble_stats(cube::YAXArray; dim="time")
    return _xmap_call(
        ensemble_stats,
        cube;
        reduced_dims=dim,
        output_axes=(Dim{:stats}(["max", "mean", "min"]),),
    )
end

ensemble_fct(cube::YAXArray; dim="time") = ensemble_stats(cube; dim=dim)

function _ensemble_mean_std_max_min_kernel(xout, xin; weights=nothing, min_members=1)
    xout .= NaN

    values, positions = _valid_numeric_values(xin)
    length(values) >= min_members || return

    if isnothing(weights)
        xout[1] = Statistics.mean(values)
        xout[2] = Statistics.std(values; corrected=false)
    else
        member_weights = Float64[weights[index] for index in positions]
        xout[1], xout[2] = _weighted_mean_std(values, member_weights)
    end

    xout[3] = Statistics.maximum(values)
    xout[4] = Statistics.minimum(values)
end

function _ensemble_mean_std_max_min_cube(cube::YAXArray; realization_dim="realization", weights=nothing, min_members=1)
    total_members = length(lookup(cube, realization_dim))
    resolved_weights = _resolve_realization_weights(weights, total_members, realization_dim)
    required_members = _required_member_count(min_members, total_members)

    return _xmap_call(
        _ensemble_mean_std_max_min_kernel,
        cube;
        reduced_dims=realization_dim,
        output_axes=(Dim{:stats}(["mean", "stdev", "max", "min"]),),
        function_kwargs=(weights=resolved_weights, min_members=required_members),
    )
end

"""
    ensemble_mean_std_max_min(data; realization_dim="realization", weights=nothing, min_members=1)

Compute xclim-style ensemble mean, standard deviation, maximum, and minimum along the realization dimension.

For a single `YAXArray`, the result is returned as a `Dataset` with variables `mean`, `stdev`, `max`, and `min`.
For a `Dataset`, each variable with a realization dimension is summarized into `<variable>_mean`, `<variable>_stdev`,
`<variable>_max`, and `<variable>_min`.
"""
function ensemble_mean_std_max_min(cube::YAXArray; realization_dim="realization", weights=nothing, min_members=1)
    summary = _ensemble_mean_std_max_min_cube(cube; realization_dim=realization_dim, weights=weights, min_members=min_members)
    stats_position = _require_axis_position(summary, :stats)

    return _dataset_from_pairs([
        :mean => _slice_axis(summary, stats_position, 1),
        :stdev => _slice_axis(summary, stats_position, 2),
        :max => _slice_axis(summary, stats_position, 3),
        :min => _slice_axis(summary, stats_position, 4),
    ])
end

function ensemble_mean_std_max_min(ds::Dataset; realization_dim="realization", weights=nothing, min_members=1)
    output_pairs = Pair{Symbol, Any}[]

    for variable_name in keys(ds.cubes)
        cube = ds[variable_name]
        _has_axis(cube, realization_dim) || continue

        summary = _ensemble_mean_std_max_min_cube(cube; realization_dim=realization_dim, weights=weights, min_members=min_members)
        stats_position = _require_axis_position(summary, :stats)

        push!(output_pairs, Symbol(string(variable_name), "_mean") => _slice_axis(summary, stats_position, 1))
        push!(output_pairs, Symbol(string(variable_name), "_stdev") => _slice_axis(summary, stats_position, 2))
        push!(output_pairs, Symbol(string(variable_name), "_max") => _slice_axis(summary, stats_position, 3))
        push!(output_pairs, Symbol(string(variable_name), "_min") => _slice_axis(summary, stats_position, 4))
    end

    isempty(output_pairs) && error("Dataset does not contain any variables with a $(realization_dim) axis.")
    return _dataset_from_pairs(output_pairs)
end

function _ensemble_percentiles_kernel(xout, xin; values, weights=nothing, min_members=1, method="linear")
    xout .= NaN

    quantile_values = Float64.(values)
    data_values, positions = _valid_numeric_values(xin)
    length(data_values) >= min_members || return

    if isnothing(weights)
        alpha, beta = _ENSEMBLE_QUANTILE_PARAMS[method]
        xout .= [Statistics.quantile(data_values, value / 100; alpha=alpha, beta=beta) for value in quantile_values]
        return
    end

    method == "linear" || error("Weighted ensemble_percentiles currently only supports the linear method.")
    member_weights = Float64[weights[index] for index in positions]
    xout .= _weighted_quantiles_1d(data_values, member_weights, quantile_values ./ 100; method=method)
end

function _ensemble_percentiles_cube(cube::YAXArray; values=[10, 50, 90], realization_dim="realization", weights=nothing, min_members=1, method="linear")
    total_members = length(lookup(cube, realization_dim))
    resolved_weights = _resolve_realization_weights(weights, total_members, realization_dim)
    required_members = _required_member_count(min_members, total_members)
    normalized_method = _normalize_quantile_method(method)
    percentile_values = Float64.(collect(values))

    return _xmap_call(
        _ensemble_percentiles_kernel,
        cube;
        reduced_dims=realization_dim,
        output_axes=(Dim{:percentiles}(percentile_values),),
        function_kwargs=(values=percentile_values, weights=resolved_weights, min_members=required_members, method=normalized_method),
    )
end

"""
    ensemble_percentiles(data; values=[10, 50, 90], realization_dim="realization", weights=nothing, min_members=1, split=true, method="linear")

Compute xclim-style ensemble percentiles along the realization dimension.

For a single `YAXArray`, `split=false` returns a single cube with a `percentiles` axis. With `split=true`, the result is
returned as a `Dataset` with variables such as `p10`, `p50`, and `p90`. For a `Dataset`, variables are either kept under
their original names with a `percentiles` axis (`split=false`) or split into `<variable>_pXX` variables (`split=true`).
"""
function ensemble_percentiles(cube::YAXArray; values=[10, 50, 90], realization_dim="realization", weights=nothing, min_members=1, split::Bool=true, method="linear")
    percentiles = _ensemble_percentiles_cube(cube; values=values, realization_dim=realization_dim, weights=weights, min_members=min_members, method=method)
    split || return percentiles

    percentile_position = _require_axis_position(percentiles, :percentiles)
    output_pairs = Pair{Symbol, Any}[]

    for (index, value) in enumerate(Float64.(collect(values)))
        push!(output_pairs, Symbol("p", _percentile_suffix(value)) => _slice_axis(percentiles, percentile_position, index))
    end

    return _dataset_from_pairs(output_pairs)
end

function ensemble_percentiles(ds::Dataset; values=[10, 50, 90], realization_dim="realization", weights=nothing, min_members=1, split::Bool=true, method="linear")
    output_pairs = Pair{Symbol, Any}[]

    for variable_name in keys(ds.cubes)
        cube = ds[variable_name]
        _has_axis(cube, realization_dim) || continue

        percentile_cube = _ensemble_percentiles_cube(cube; values=values, realization_dim=realization_dim, weights=weights, min_members=min_members, method=method)

        if split
            percentile_position = _require_axis_position(percentile_cube, :percentiles)
            for (index, value) in enumerate(Float64.(collect(values)))
                push!(output_pairs, Symbol(string(variable_name), "_p", _percentile_suffix(value)) => _slice_axis(percentile_cube, percentile_position, index))
            end
        else
            push!(output_pairs, variable_name => percentile_cube)
        end
    end

    isempty(output_pairs) && error("Dataset does not contain any variables with a $(realization_dim) axis.")
    return _dataset_from_pairs(output_pairs)
end

function _criteria_labels(cube::YAXArray, other_axes::Vector{Int}, prefix)
    if isempty(other_axes)
        return [isnothing(prefix) ? "value" : string(prefix)]
    end

    axis_names = _axis_names(cube)
    coordinate_sets = [collect(lookup(cube, axis_names[index])) for index in other_axes]
    labels = String[]

    for coordinates in Iterators.product(coordinate_sets...)
        parts = String[]
        isnothing(prefix) || push!(parts, string(prefix))
        for (axis_index, coordinate) in zip(other_axes, coordinates)
            push!(parts, string(axis_names[axis_index], "=", coordinate))
        end
        push!(labels, join(parts, "|"))
    end

    return labels
end

function _flatten_criteria(cube::YAXArray; realization_dim="realization", prefix=nothing)
    realization_position = _require_axis_position(cube, realization_dim)
    axis_indices = collect(eachindex(cube.axes))
    other_axes = [index for index in axis_indices if index != realization_position]
    order = (realization_position, other_axes...)

    data = map(_float_or_nan, permutedims(Array(cube), order))
    matrix = reshape(data, size(data, 1), :)
    labels = _criteria_labels(cube, other_axes, prefix)

    keep_mask = [any(value -> !isnan(value), matrix[:, column]) for column in axes(matrix, 2)]
    return matrix[:, keep_mask], labels[keep_mask]
end

"""
    make_criteria(data; realization_dim="realization")

Flatten all non-realization dimensions into a `criteria` axis, following xclim's ensemble-reduction preprocessing.
For `Dataset` inputs, each variable is flattened independently and concatenated along the output criteria axis.
Criteria columns that are entirely NaN are removed.
"""
function make_criteria(cube::YAXArray; realization_dim="realization")
    matrix, labels = _flatten_criteria(cube; realization_dim=realization_dim)
    realization_symbol = _normalize_xmap_dim(realization_dim)
    return YAXArray((Dim{realization_symbol}(lookup(cube, realization_symbol)), Dim{:criteria}(labels)), matrix)
end

function make_criteria(ds::Dataset; realization_dim="realization")
    matrix_blocks = Matrix{Float64}[]
    labels = String[]
    realization_lookup = nothing

    for variable_name in keys(ds.cubes)
        cube = ds[variable_name]
        _has_axis(cube, realization_dim) || continue

        block, block_labels = _flatten_criteria(cube; realization_dim=realization_dim, prefix=variable_name)
        isempty(block_labels) && continue

        if isnothing(realization_lookup)
            realization_lookup = lookup(cube, realization_dim)
        elseif collect(realization_lookup) != collect(lookup(cube, realization_dim))
            error("All variables passed to make_criteria must share the same realization axis.")
        end

        push!(matrix_blocks, block)
        append!(labels, block_labels)
    end

    isempty(matrix_blocks) && error("Dataset does not contain any variables with a $(realization_dim) axis.")
    matrix = hcat(matrix_blocks...)
    realization_symbol = _normalize_xmap_dim(realization_dim)
    return YAXArray((Dim{realization_symbol}(realization_lookup), Dim{:criteria}(labels)), matrix)
end

function _standardize_matrix(matrix::AbstractMatrix{<:Real})
    mean_values = vec(mean(matrix; dims=1))
    std_values = vec(std(matrix; dims=1, corrected=true))
    std_values[.!isfinite.(std_values) .| (std_values .== 0)] .= 1.0
    standardized = (Float64.(matrix) .- reshape(mean_values, 1, :)) ./ reshape(std_values, 1, :)
    return standardized, mean_values, std_values
end

function _covariance_matrix(matrix::AbstractMatrix{<:Real})
    centered = Float64.(matrix) .- mean(matrix; dims=1)
    denominator = max(size(centered, 1) - 1, 1)
    return (centered' * centered) / denominator
end

function _distance_between(x::AbstractVector, y::AbstractVector, dist_method::AbstractString; kwargs...)
    dx = Float64.(x)
    dy = Float64.(y)
    delta = dx .- dy
    normalized = lowercase(dist_method)

    if normalized == "euclidean"
        return norm(delta)
    elseif normalized == "sqeuclidean"
        return sum(abs2, delta)
    elseif normalized == "cityblock"
        return sum(abs, delta)
    elseif normalized == "chebyshev"
        return maximum(abs, delta)
    elseif normalized == "cosine"
        denominator = norm(dx) * norm(dy)
        return iszero(denominator) ? 0.0 : 1 - dot(dx, dy) / denominator
    elseif normalized == "seuclidean"
        variances = Float64.(get(kwargs, :V, ones(length(delta))))
        return sqrt(sum((delta .^ 2) ./ variances))
    elseif normalized == "mahalanobis"
        inverse_covariance = get(kwargs, :VI, nothing)
        isnothing(inverse_covariance) && error("VI must be provided for mahalanobis distance.")

        if inverse_covariance isa AbstractVector
            return sqrt(sum((delta .^ 2) .* Float64.(inverse_covariance)))
        end

        matrix = Float64.(Matrix(inverse_covariance))
        return sqrt(dot(delta, matrix, delta))
    end

    error("Unsupported distance method $(dist_method).")
end

"""
    kkz_reduce_ensemble(data, num_select; realization_dim="realization", criteria_dim="criteria", dist_method="euclidean", standardize=true, kwargs...)

Select a representative subset of realizations using the deterministic KKZ algorithm.
If `data` is not already a realization-by-criteria matrix, `make_criteria` is applied first.
"""
function kkz_reduce_ensemble(data::YAXArray, num_select::Integer; realization_dim="realization", criteria_dim="criteria", dist_method="euclidean", standardize::Bool=true, kwargs...)
    criteria_cube = if _has_axis(data, criteria_dim) && ndims(data) == 2
        data
    else
        make_criteria(data; realization_dim=realization_dim)
    end

    matrix = Float64.(Array(criteria_cube))
    any(isnan, matrix) && error("kkz_reduce_ensemble requires criteria without NaN values.")

    total_members = size(matrix, 1)
    1 <= num_select <= total_members || error("num_select must be between 1 and the number of realizations.")

    processed = copy(matrix)
    metric_kwargs = Dict{Symbol, Any}(kwargs)

    if standardize
        processed, _, _ = _standardize_matrix(processed)
    end

    normalized_metric = lowercase(String(dist_method))
    if normalized_metric == "seuclidean" && !haskey(metric_kwargs, :V)
        metric_kwargs[:V] = if standardize
            ones(size(processed, 2))
        else
            variances = vec(var(processed; dims=1, corrected=true))
            variances[.!isfinite.(variances) .| (variances .== 0)] .= 1.0
            variances
        end
    elseif normalized_metric == "mahalanobis" && !haskey(metric_kwargs, :VI)
        metric_kwargs[:VI] = pinv(_covariance_matrix(processed))
    end

    selected = Int[]
    remaining = collect(1:total_members)
    centroid = vec(mean(processed; dims=1))

    centroid_distances = [
        _distance_between(centroid, view(processed, candidate, :), normalized_metric; metric_kwargs...)
        for candidate in remaining
    ]
    first_position = argmin(centroid_distances)
    push!(selected, remaining[first_position])
    deleteat!(remaining, first_position)

    while length(selected) < num_select
        separation = Float64[]
        for candidate in remaining
            candidate_vector = view(processed, candidate, :)
            push!(
                separation,
                minimum(_distance_between(view(processed, chosen, :), candidate_vector, normalized_metric; metric_kwargs...) for chosen in selected),
            )
        end

        chosen_position = argmax(separation)
        push!(selected, remaining[chosen_position])
        deleteat!(remaining, chosen_position)
    end

    return selected
end

function kkz_reduce_ensemble(ds::Dataset, num_select::Integer; realization_dim="realization", kwargs...)
    return kkz_reduce_ensemble(make_criteria(ds; realization_dim=realization_dim), num_select; realization_dim=realization_dim, kwargs...)
end

function _outer_axes(cube::YAXArray, excluded_dims)
    axis_names = _axis_names(cube)
    outer_positions = [index for index in eachindex(axis_names) if !(axis_names[index] in excluded_dims)]
    return Tuple(cube.axes[index] for index in outer_positions), [axis_names[index] for index in outer_positions]
end

function _flatten_delta_cube(cube::YAXArray, outer_dims::Vector{Symbol}, realization_dim::Symbol)
    axis_names = _axis_names(cube)
    order = Int[]
    append!(order, (_require_axis_position(cube, dim) for dim in outer_dims))
    has_realization = realization_dim in axis_names
    has_realization && push!(order, _require_axis_position(cube, realization_dim))

    data = isempty(order) ? map(_float_or_nan, Array(cube)) : map(_float_or_nan, permutedims(Array(cube), Tuple(order)))
    outer_sizes = size(data)[1:length(outer_dims)]

    if has_realization
        flattened = reshape(data, :, size(data, length(outer_dims) + 1))
    else
        flattened = reshape(data, :, 1)
    end

    return flattened, outer_sizes
end

function _flatten_timeseries_cube(cube::YAXArray, outer_dims::Vector{Symbol}, realization_dim::Symbol, time_dim::Symbol)
    axis_names = _axis_names(cube)
    order = Int[]
    append!(order, (_require_axis_position(cube, dim) for dim in outer_dims))
    has_realization = realization_dim in axis_names
    has_realization && push!(order, _require_axis_position(cube, realization_dim))
    push!(order, _require_axis_position(cube, time_dim))

    data = map(_float_or_nan, permutedims(Array(cube), Tuple(order)))
    outer_sizes = size(data)[1:length(outer_dims)]

    if has_realization
        flattened = reshape(data, :, size(data, length(outer_dims) + 1), size(data, length(outer_dims) + 2))
    else
        flattened = reshape(data, :, 1, size(data, length(outer_dims) + 1))
    end

    return flattened, outer_sizes
end

function _series_is_valid(series::AbstractVector{<:Real}, invalid)
    valid_count = count(!isnan, series)
    if isnothing(invalid)
        return valid_count == length(series)
    end

    required = if invalid isa Integer
        Int(invalid)
    elseif hasproperty(invalid, :n)
        Int(getproperty(invalid, :n))
    else
        error("invalid must be nothing, an integer minimum count, or an object exposing an n field.")
    end

    return valid_count >= required
end

function _detrend_series(series::AbstractVector{<:Real}, degree::Integer)
    valid = findall(!isnan, series)
    length(valid) > degree || return fill(NaN, length(series))

    x = Float64.(valid)
    y = Float64.(series[valid])
    design = hcat((x .^ power for power in 0:degree)...)
    coefficients = design \ y
    fitted = design * coefficients
    residuals = fill(NaN, length(series))
    residuals[valid] = y .- fitted
    return residuals
end

function _pvalue_ttest(future_series::AbstractVector{<:Real}, reference_series::AbstractVector{<:Real})
    future_values = Float64.(future_series[.!isnan.(future_series)])
    reference_values = Float64.(reference_series[.!isnan.(reference_series)])
    length(future_values) >= 2 && !isempty(reference_values) || return NaN

    sigma = Statistics.std(future_values; corrected=true)
    iszero(sigma) && return NaN
    t_stat = (Statistics.mean(future_values) - Statistics.mean(reference_values)) / (sigma / sqrt(length(future_values)))
    return 2 * ccdf(TDist(length(future_values) - 1), abs(t_stat))
end

function _pvalue_welch_ttest(future_series::AbstractVector{<:Real}, reference_series::AbstractVector{<:Real})
    future_values = Float64.(future_series[.!isnan.(future_series)])
    reference_values = Float64.(reference_series[.!isnan.(reference_series)])
    length(future_values) >= 2 && length(reference_values) >= 2 || return NaN

    future_variance = Statistics.var(future_values; corrected=true)
    reference_variance = Statistics.var(reference_values; corrected=true)
    scale = future_variance / length(future_values) + reference_variance / length(reference_values)
    iszero(scale) && return NaN

    t_stat = (Statistics.mean(future_values) - Statistics.mean(reference_values)) / sqrt(scale)
    numerator = scale^2
    denominator = future_variance^2 / (length(future_values)^2 * (length(future_values) - 1)) + reference_variance^2 / (length(reference_values)^2 * (length(reference_values) - 1))
    iszero(denominator) && return NaN
    degrees_of_freedom = numerator / denominator
    return 2 * ccdf(TDist(degrees_of_freedom), abs(t_stat))
end

function _average_ranks(values::AbstractVector{<:Real})
    order = sortperm(values)
    sorted_values = values[order]
    ranks = similar(Float64.(values))
    index = 1

    while index <= length(values)
        upper = index
        while upper < length(values) && sorted_values[upper + 1] == sorted_values[index]
            upper += 1
        end

        average_rank = (index + upper) / 2
        for position in index:upper
            ranks[order[position]] = average_rank
        end
        index = upper + 1
    end

    return ranks
end

function _pvalue_mannwhitney(future_series::AbstractVector{<:Real}, reference_series::AbstractVector{<:Real})
    future_values = Float64.(future_series[.!isnan.(future_series)])
    reference_values = Float64.(reference_series[.!isnan.(reference_series)])
    length(future_values) >= 1 && length(reference_values) >= 1 || return NaN

    combined = vcat(future_values, reference_values)
    ranks = _average_ranks(combined)
    n_future = length(future_values)
    n_reference = length(reference_values)
    rank_sum = sum(ranks[1:n_future])
    u_statistic = rank_sum - n_future * (n_future + 1) / 2
    mean_u = n_future * n_reference / 2
    variance_u = n_future * n_reference * (n_future + n_reference + 1) / 12
    iszero(variance_u) && return NaN

    z = (u_statistic - mean_u) / sqrt(variance_u)
    return 2 * ccdf(Normal(), abs(z))
end

function _pvalue_brownforsythe(future_series::AbstractVector{<:Real}, reference_series::AbstractVector{<:Real})
    future_values = Float64.(future_series[.!isnan.(future_series)])
    reference_values = Float64.(reference_series[.!isnan.(reference_series)])
    length(future_values) >= 2 && length(reference_values) >= 2 || return NaN

    future_deviation = abs.(future_values .- Statistics.median(future_values))
    reference_deviation = abs.(reference_values .- Statistics.median(reference_values))
    pooled = vcat(future_deviation, reference_deviation)
    overall_mean = Statistics.mean(pooled)
    group_means = (Statistics.mean(future_deviation), Statistics.mean(reference_deviation))

    numerator = length(future_deviation) * (group_means[1] - overall_mean)^2 + length(reference_deviation) * (group_means[2] - overall_mean)^2
    denominator = sum((future_deviation .- group_means[1]) .^ 2) + sum((reference_deviation .- group_means[2]) .^ 2)
    denominator /= (length(pooled) - 2)
    iszero(denominator) && return NaN

    statistic = numerator / denominator
    return ccdf(FDist(1, length(pooled) - 2), statistic)
end

function _pvalue_for_test(future_series::AbstractVector{<:Real}, reference_series::AbstractVector{<:Real}, test)
    normalized_test = lowercase(String(test))

    if normalized_test == "ttest"
        return _pvalue_ttest(future_series, reference_series)
    elseif normalized_test == "welch-ttest"
        return _pvalue_welch_ttest(future_series, reference_series)
    elseif normalized_test == "mannwhitney-utest"
        return _pvalue_mannwhitney(future_series, reference_series)
    elseif normalized_test == "brownforsythe-test"
        return _pvalue_brownforsythe(future_series, reference_series)
    end

    error("Unsupported robustness test $(test).")
end

function _ipcc_ar6_gamma(reference_series::AbstractVector{<:Real}; ref_pi=nothing)
    if isnothing(ref_pi)
        detrended = _detrend_series(reference_series, 1)
        valid = detrended[.!isnan.(detrended)]
        isempty(valid) && return NaN
        return sqrt(2 / 20) * 1.645 * Statistics.std(valid; corrected=false)
    end

    detrended = _detrend_series(ref_pi, 2)
    valid = detrended[.!isnan.(detrended)]
    length(valid) >= 20 || return NaN
    chunk_count = div(length(valid), 20)
    chunk_means = [Statistics.mean(valid[(20 * (index - 1) + 1):(20 * index)]) for index in 1:chunk_count]
    isempty(chunk_means) && return NaN
    return sqrt(2) * 1.645 * Statistics.std(chunk_means; corrected=false)
end

function _comparison_result(future_series::AbstractVector{<:Real}, reference_series::AbstractVector{<:Real}, delta::Float64, test, kwargs)
    isnothing(test) && return true, NaN

    normalized_test = lowercase(String(test))
    if normalized_test == "threshold"
        abs_threshold = get(kwargs, :abs_thresh, nothing)
        rel_threshold = get(kwargs, :rel_thresh, nothing)

        if !isnothing(abs_threshold) && isnothing(rel_threshold)
            return abs(delta) > abs_threshold, NaN
        elseif !isnothing(rel_threshold) && isnothing(abs_threshold)
            reference_mean = Statistics.mean(Float64.(reference_series[.!isnan.(reference_series)]))
            iszero(reference_mean) && return false, NaN
            return abs(delta / reference_mean) > rel_threshold, NaN
        end

        error("One and only one of abs_thresh or rel_thresh must be supplied when test='threshold'.")
    elseif normalized_test == "ipcc-ar6-c"
        gamma = _ipcc_ar6_gamma(reference_series; ref_pi=get(kwargs, :ref_pi, nothing))
        return isfinite(gamma) && abs(delta) > gamma, NaN
    end

    p_value = _pvalue_for_test(future_series, reference_series, normalized_test)
    return isfinite(p_value) && p_value < get(kwargs, :p_change, 0.05), p_value
end

function _output_cube(axes, data, properties=Dict{String, Any}())
    return YAXArray(axes, data, properties)
end

"""
    robustness_fractions(fut, ref=nothing; test=nothing, weights=nothing, invalid=nothing, realization_dim="realization", time_dim="time", kwargs...)

Compute xclim-style robustness fractions describing the existence and sign of change across an ensemble.

`fut` should be a future-period cube with `realization` and `time` dimensions. If `ref` is omitted, `fut` is treated as
delta values along the realization dimension and only `test=nothing` or `test="threshold"` are allowed.

The result is a `Dataset` with variables `changed`, `positive`, `changed_positive`, `negative`, `changed_negative`,
`agree`, and `valid`. When a significance test produces p-values, a `pvals` variable is also returned.
"""
function robustness_fractions(fut::YAXArray, ref::Union{Nothing, YAXArray}=nothing; test=nothing, weights=nothing, invalid=nothing, realization_dim="realization", time_dim="time", kwargs...)
    realization_symbol = _normalize_xmap_dim(realization_dim)
    time_symbol = _normalize_xmap_dim(time_dim)
    fut_axis_names = _axis_names(fut)
    outer_dims = [dim for dim in fut_axis_names if dim != realization_symbol && dim != time_symbol]
    outer_axes, _ = _outer_axes(fut, Set([realization_symbol, time_symbol]))

    if isnothing(ref)
        lower_test = isnothing(test) ? nothing : lowercase(String(test))
        lower_test in (nothing, "threshold") || error("When ref is omitted, test must be nothing or 'threshold'.")

        fut_values, outer_sizes = _flatten_delta_cube(fut, outer_dims, realization_symbol)
        total_members = size(fut_values, 2)
        resolved_weights = _resolve_realization_weights(weights, total_members, realization_symbol)
        member_weights = isnothing(resolved_weights) ? ones(Float64, total_members) : resolved_weights

        output_shape = isempty(outer_sizes) ? () : Tuple(outer_sizes)
        output_length = size(fut_values, 1)
        changed_fraction = fill(0.0, output_length)
        positive_fraction = fill(0.0, output_length)
        changed_positive_fraction = fill(0.0, output_length)
        negative_fraction = fill(0.0, output_length)
        changed_negative_fraction = fill(0.0, output_length)
        agree_fraction = fill(0.0, output_length)
        valid_fraction = fill(0.0, output_length)

        for index in 1:size(fut_values, 1)
            delta = view(fut_values, index, :)
            valid = .!isnan.(delta)
            valid_weight = sum(member_weights[valid])
            valid_fraction[index] = valid_weight / total_members

            if iszero(valid_weight)
                continue
            end

            changed = if isnothing(test)
                trues(total_members)
            else
                abs_threshold = get(kwargs, :abs_thresh, nothing)
                rel_threshold = get(kwargs, :rel_thresh, nothing)
                if !isnothing(abs_threshold) && isnothing(rel_threshold)
                    abs.(delta) .> abs_threshold
                elseif !isnothing(rel_threshold) && isnothing(abs_threshold)
                    abs.(delta) .> rel_threshold
                else
                    error("One and only one of abs_thresh or rel_thresh must be supplied when test='threshold'.")
                end
            end

            positive = delta .> 0
            negative = delta .< 0

            changed_fraction[index] = sum(member_weights[valid .& changed]) / valid_weight
            positive_fraction[index] = sum(member_weights[valid .& positive]) / valid_weight
            negative_fraction[index] = sum(member_weights[valid .& negative]) / valid_weight
            changed_positive_fraction[index] = sum(member_weights[valid .& changed .& positive]) / valid_weight
            changed_negative_fraction[index] = sum(member_weights[valid .& changed .& negative]) / valid_weight
            agree_fraction[index] = max(positive_fraction[index], negative_fraction[index], 1 - positive_fraction[index] - negative_fraction[index])
        end

        output_axes_tuple = outer_axes
        output_pairs = Pair{Symbol, Any}[
            :changed => _output_cube(output_axes_tuple, _reshape_output(changed_fraction, output_shape), Dict("test" => string(test))),
            :positive => _output_cube(output_axes_tuple, _reshape_output(positive_fraction, output_shape)),
            :changed_positive => _output_cube(output_axes_tuple, _reshape_output(changed_positive_fraction, output_shape), Dict("test" => string(test))),
            :negative => _output_cube(output_axes_tuple, _reshape_output(negative_fraction, output_shape)),
            :changed_negative => _output_cube(output_axes_tuple, _reshape_output(changed_negative_fraction, output_shape), Dict("test" => string(test))),
            :agree => _output_cube(output_axes_tuple, _reshape_output(agree_fraction, output_shape)),
            :valid => _output_cube(output_axes_tuple, _reshape_output(valid_fraction, output_shape)),
        ]
        return _dataset_from_pairs(output_pairs)
    end

    fut_values, outer_sizes = _flatten_timeseries_cube(fut, outer_dims, realization_symbol, time_symbol)
    ref_values, _ = _flatten_timeseries_cube(ref, outer_dims, realization_symbol, time_symbol)
    total_members = size(fut_values, 2)
    resolved_weights = _resolve_realization_weights(weights, total_members, realization_symbol)
    member_weights = isnothing(resolved_weights) ? ones(Float64, total_members) : resolved_weights
    ref_member_count = size(ref_values, 2)

    output_shape = isempty(outer_sizes) ? () : Tuple(outer_sizes)
    output_length = size(fut_values, 1)
    changed_fraction = fill(0.0, output_length)
    positive_fraction = fill(0.0, output_length)
    changed_positive_fraction = fill(0.0, output_length)
    negative_fraction = fill(0.0, output_length)
    changed_negative_fraction = fill(0.0, output_length)
    agree_fraction = fill(0.0, output_length)
    valid_fraction = fill(0.0, output_length)
    pvals = fill(NaN, output_length, total_members)
    has_pvals = false

    for index in 1:size(fut_values, 1)
        delta = fill(NaN, total_members)
        changed = falses(total_members)
        positive = falses(total_members)
        negative = falses(total_members)
        valid = falses(total_members)

        for member in 1:total_members
            future_series = vec(fut_values[index, member, :])
            reference_series = if ref_member_count == 1
                vec(ref_values[index, 1, :])
            else
                vec(ref_values[index, member, :])
            end

            valid[member] = _series_is_valid(future_series, invalid) && _series_is_valid(reference_series, invalid)
            valid[member] || continue

            delta[member] = Statistics.mean(Float64.(future_series[.!isnan.(future_series)])) - Statistics.mean(Float64.(reference_series[.!isnan.(reference_series)]))
            changed[member], p_value = _comparison_result(future_series, reference_series, delta[member], test, kwargs)
            positive[member] = delta[member] > 0
            negative[member] = delta[member] < 0

            if isfinite(p_value)
                pvals[index, member] = p_value
                has_pvals = true
            end
        end

        valid_weight = sum(member_weights[valid])
        valid_fraction[index] = valid_weight / total_members

        if iszero(valid_weight)
            continue
        end

        changed_fraction[index] = sum(member_weights[valid .& changed]) / valid_weight
        positive_fraction[index] = sum(member_weights[valid .& positive]) / valid_weight
        negative_fraction[index] = sum(member_weights[valid .& negative]) / valid_weight
        changed_positive_fraction[index] = sum(member_weights[valid .& changed .& positive]) / valid_weight
        changed_negative_fraction[index] = sum(member_weights[valid .& changed .& negative]) / valid_weight
        agree_fraction[index] = max(positive_fraction[index], negative_fraction[index], 1 - positive_fraction[index] - negative_fraction[index])
    end

    output_axes_tuple = outer_axes
    output_pairs = Pair{Symbol, Any}[
        :changed => _output_cube(output_axes_tuple, _reshape_output(changed_fraction, output_shape), Dict("test" => string(test))),
        :positive => _output_cube(output_axes_tuple, _reshape_output(positive_fraction, output_shape)),
        :changed_positive => _output_cube(output_axes_tuple, _reshape_output(changed_positive_fraction, output_shape), Dict("test" => string(test))),
        :negative => _output_cube(output_axes_tuple, _reshape_output(negative_fraction, output_shape)),
        :changed_negative => _output_cube(output_axes_tuple, _reshape_output(changed_negative_fraction, output_shape), Dict("test" => string(test))),
        :agree => _output_cube(output_axes_tuple, _reshape_output(agree_fraction, output_shape)),
        :valid => _output_cube(output_axes_tuple, _reshape_output(valid_fraction, output_shape)),
    ]

    if has_pvals
        pval_shape = isempty(output_shape) ? (total_members,) : (output_shape..., total_members)
        pval_axes = isempty(output_axes_tuple) ? (Dim{realization_symbol}(1:total_members),) : (output_axes_tuple..., Dim{realization_symbol}(1:total_members))
        push!(output_pairs, :pvals => _output_cube(pval_axes, reshape(copy(pvals), pval_shape)))
    end

    return _dataset_from_pairs(output_pairs)
end

function _compare_fraction(value, op, threshold)
    if isnothing(op) || isempty(String(op))
        return true
    end

    isnan(value) && return false

    if op == ">="
        return value >= threshold
    elseif op == ">"
        return value > threshold
    elseif op == "<="
        return value <= threshold
    elseif op == "<"
        return value < threshold
    elseif op == "=="
        return value == threshold
    end

    error("Unsupported comparison operator $(op).")
end

"""
    robustness_categories(changed_or_fractions, agree=nothing, valid=nothing; categories=nothing, ops=nothing, thresholds=nothing)

Create an IPCC-style categorical robustness map from robustness fractions.
The default thresholds correspond to xclim's AR6-aligned defaults.
"""
function robustness_categories(changed_or_fractions, agree=nothing, valid=nothing; categories=nothing, ops=nothing, thresholds=nothing)
    if categories === nothing
        categories = [
            "Robust signal",
            "No change or no signal",
            "Conflicting signal",
        ]
    end

    if ops === nothing
        ops = [(">=", ">="), ("<", nothing), (">=", "<")]
    end

    if thresholds === nothing
        thresholds = [(0.66, 0.8), (0.66, nothing), (0.66, 0.8)]
    end

    changed = if changed_or_fractions isa Dataset
        ds = changed_or_fractions
        isnothing(agree) && (agree = ds[:agree])
        isnothing(valid) && (valid = ds[:valid])
        ds[:changed]
    else
        changed_or_fractions
    end

    changed_data = map(_float_or_nan, Array(changed))
    agree_data = isnothing(agree) ? fill(NaN, size(changed_data)) : map(_float_or_nan, Array(agree))
    output = fill(99, size(changed_data))

    for index in eachindex(changed_data)
        for category_index in eachindex(categories)
            change_rule, agree_rule = ops[category_index]
            change_threshold, agree_threshold = thresholds[category_index]
            if _compare_fraction(changed_data[index], change_rule, change_threshold) && _compare_fraction(agree_data[index], agree_rule, agree_threshold)
                output[index] = category_index
                break
            end
        end
    end

    if !isnothing(valid)
        valid_data = map(_float_or_nan, Array(valid))
        output[valid_data .<= 0] .= 99
    end

    properties = Dict{String, Any}(
        "flag_values" => collect(1:length(categories)),
        "_FillValue" => 99,
        "flag_descriptions" => categories,
        "flag_meanings" => join(map(category -> replace(lowercase(category), " " => "_"), categories), " "),
    )
    return _output_cube(changed.axes, output, properties)
end

function _cdf_squared_area(x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real})
    values1 = sort(Float64.(x1[.!isnan.(x1)]))
    values2 = sort(Float64.(x2[.!isnan.(x2)]))
    (isempty(values1) || isempty(values2)) && return NaN

    grid = sort(unique(vcat(values1, values2)))
    length(grid) <= 1 && return 0.0
    area = 0.0

    for index in 1:(length(grid) - 1)
        left = grid[index]
        right = grid[index + 1]
        cdf1 = searchsortedlast(values1, left) / length(values1)
        cdf2 = searchsortedlast(values2, left) / length(values2)
        area += (cdf1 - cdf2)^2 * (right - left)
    end

    return area
end

function _robustness_coefficient_arrays(reference_series::AbstractVector{<:Real}, future_matrix::AbstractMatrix{<:Real})
    reference_values = Float64.(reference_series[.!isnan.(reference_series)])
    future_values = map(_float_or_nan, future_matrix)
    valid_future = future_values[:, vec(any(.!isnan.(future_values); dims=1))]
    isempty(reference_values) && return NaN
    size(valid_future, 2) == 0 && return NaN

    pooled_future = Float64[valid_future[index] for index in eachindex(valid_future) if !isnan(valid_future[index])]
    model_mean_projection = vec(mean(valid_future; dims=1))
    model_mean_projection = model_mean_projection[.!isnan.(model_mean_projection)]
    isempty(pooled_future) && return NaN
    isempty(model_mean_projection) && return NaN

    area_future = _cdf_squared_area(pooled_future, model_mean_projection)
    area_reference = _cdf_squared_area(reference_values, model_mean_projection)
    iszero(area_reference) && return NaN
    return 1 - area_future / area_reference
end

"""
    robustness_coefficient(fut, ref; realization_dim="realization", time_dim="time")

Compute the Knutti-Sedlacek robustness coefficient for an ensemble signal.
"""
function robustness_coefficient(fut::YAXArray, ref::YAXArray; realization_dim="realization", time_dim="time")
    realization_symbol = _normalize_xmap_dim(realization_dim)
    time_symbol = _normalize_xmap_dim(time_dim)
    outer_dims = [dim for dim in _axis_names(fut) if dim != realization_symbol && dim != time_symbol]
    outer_axes, _ = _outer_axes(fut, Set([realization_symbol, time_symbol]))
    future_values, outer_sizes = _flatten_timeseries_cube(fut, outer_dims, realization_symbol, time_symbol)
    reference_values, _ = _flatten_timeseries_cube(ref, outer_dims, realization_symbol, time_symbol)
    ref_member_count = size(reference_values, 2)
    output_shape = isempty(outer_sizes) ? () : Tuple(outer_sizes)
    coefficients = fill(NaN, size(future_values, 1))

    for index in 1:size(future_values, 1)
        reference_series = if ref_member_count == 1
            vec(reference_values[index, 1, :])
        else
            vec(mean(reference_values[index, :, :]; dims=1))
        end
        coefficients[index] = _robustness_coefficient_arrays(reference_series, future_values[index, :, :])
    end

    return _output_cube(outer_axes, _reshape_output(coefficients, output_shape), Dict("name" => "R"))
end

function robustness_coefficient(fut::Dataset, ref::Dataset; realization_dim="realization", time_dim="time")
    output_pairs = Pair{Symbol, Any}[]
    for variable_name in intersect(Set(keys(fut.cubes)), Set(keys(ref.cubes)))
        future_cube = fut[variable_name]
        reference_cube = ref[variable_name]
        _has_axis(future_cube, realization_dim) || continue
        push!(output_pairs, variable_name => robustness_coefficient(future_cube, reference_cube; realization_dim=realization_dim, time_dim=time_dim))
    end

    isempty(output_pairs) && error("No shared dataset variables with a $(realization_dim) axis were found.")
    return _dataset_from_pairs(output_pairs)
end
