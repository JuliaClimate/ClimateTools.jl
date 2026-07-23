function rlevels_cube(xout, xin; threshold=nothing, rlevels = [1, 2, 5, 10, 25, 50, 100, 1000], minimalvalue=1.0)
    values = _finite_real_values(xin)
    if isempty(values)
        xout .= missing
        return
    end

    threshold_value = _gp_threshold(values; threshold=threshold, threshold_quantile=0.92, minimalvalue=minimalvalue)
    if ismissing(threshold_value)
        xout .= missing
        return
    end

    exceedances = values[values .> threshold_value] .- threshold_value
    if isempty(exceedances)
        xout .= missing
        return
    end

    try
        model = Extremes.gpfit(exceedances)
        nobs = size(xin, 1)
        nobsperblock = 365
        r_h = returnlevel.(model, threshold_value, nobs, nobsperblock, rlevels)
        xout .= [x.value[1] for x in r_h]
    catch err
        if err isa MethodError || err isa UndefKeywordError
            rethrow(err)
        end
        xout .= missing
    end

end


"""
    rlevels_cube(ds::YAXArray; rlevels = [1, 2, 5, 10, 25, 50, 100, 1000], threshold=nothing, minimalvalue=1.0, latsouth_pixel=0, latnorth_pixel=0, lonwest_pixel=0, loneast_pixel=0)

Compute the return levels for a dataset `ds` along the time dimension.

# Arguments
- `ds`: The input dataset.
- `rlevels`: The return levels to compute. Default is [1, 2, 5, 10, 25, 50, 100, 1000].
- `threshold`: The threshold to use for the Generalized Pareto fit. If not provided, the 0.92 quantile is used.
- `minimalvalue`: The minimal value to consider in the dataset. Default is 1.0.
- `latsouth_pixel`: The south pixel to consider. Default is 0.
- `latnorth_pixel`: The north pixel to consider. Default is 0.
- `lonwest_pixel`: The west pixel to consider. Default is 0.
- `loneast_pixel`: The east pixel to consider. Default is 0.

# Returns
The return levels of the input dataset `ds` along the time dimension.

"""
function rlevels_cube(ds::YAXArray; threshold=nothing, rlevels = [1, 2, 5, 10, 25, 50, 100, 1000], minimalvalue=1.0, latsouth_pixel=0, latnorth_pixel=0, lonwest_pixel=0, loneast_pixel=0, lonname="longitude", latname="latitude")
    return _xmap_call(
        rlevels_cube,
        ds;
        reduced_dims=:time,
        output_axes=(Dim{:rlevels}(rlevels),),
        function_kwargs=(threshold=threshold, rlevels=rlevels, minimalvalue=minimalvalue),
        outtype=Union{Missing, Float64},
    )

end

# function moving_fct(cube::YAXArray; fct::Function=mean, latsouth_pixel=1, latnorth_pixel=1, lonwest_pixel=1, loneast_pixel=1)
    
#     indims = InDims(MovingWindow("latitude", latsouth_pixel,latnorth_pixel),
#         MovingWindow("longitude", lonwest_pixel,loneast_pixel))
#     outdims=OutDims()
#     mapCube(moving_fct, cube; fct=fct, indims=indims, outdims=outdims)
# end


const _EXTREME_MODEL_OUTPUT = Union{Missing, Extremes.MaximumLikelihoodAbstractExtremeValueModel}

function _extreme_axis_names(cube::YAXArray)
    return Tuple(Symbol(name(axis)) for axis in cube.axes)
end

function _extreme_output_layout(cube::YAXArray, reduced_dims)
    normalized_dims = _normalize_xmap_dims(reduced_dims)
    axis_names = _extreme_axis_names(cube)

    reduced_positions = map(normalized_dims) do dim
        pos = findfirst(==(dim), axis_names)
        pos === nothing && error("Reduced dimension $(dim) not found in cube axes $(axis_names).")
        pos
    end

    keep_positions = Tuple(index for index in eachindex(cube.axes) if !(index in reduced_positions))
    isempty(keep_positions) && error("At least one non-reduced dimension must remain when fitting extreme value models on a cube.")

    output_axes = Tuple(cube.axes[index] for index in keep_positions)
    output_size = Tuple(size(cube, index) for index in keep_positions)
    return normalized_dims, output_axes, keep_positions, output_size
end

function _extreme_slice(cube::YAXArray, keep_positions::Tuple, idx::CartesianIndex)
    selectors = Any[Colon() for _ in 1:ndims(cube)]

    for (offset, pos) in enumerate(keep_positions)
        selectors[pos] = idx[offset]
    end

    return view(cube, selectors...)
end

function _finite_real_values(data)
    values = Float64[]

    for value in data
        ismissing(value) && continue
        value isa Real || continue

        numeric = Float64(value)
        isfinite(numeric) || continue
        push!(values, numeric)
    end

    return values
end

function _fit_extreme_model_cube(fitpoint::Function, ds::YAXArray; reduced_dims=:time)
    _, output_axes, keep_positions, output_size = _extreme_output_layout(ds, reduced_dims)
    models = Array{_EXTREME_MODEL_OUTPUT}(undef, output_size...)
    fill!(models, missing)

    for idx in CartesianIndices(models)
        models[idx] = fitpoint(_extreme_slice(ds, keep_positions, idx))
    end

    return YAXArray(output_axes, models)
end

function _extreme_n_observations(cube::YAXArray, reduced_dims)
    normalized_dims = _normalize_xmap_dims(reduced_dims)
    axis_names = _extreme_axis_names(cube)

    reduced_positions = map(normalized_dims) do dim
        pos = findfirst(==(dim), axis_names)
        pos === nothing && error("Reduced dimension $(dim) not found in cube axes $(axis_names).")
        pos
    end

    return Int(prod(size(cube, pos) for pos in reduced_positions))
end

function _gp_threshold(values::Vector{Float64}; threshold=nothing, threshold_quantile::Real=0.95, minimalvalue::Real=1.0)
    if isnothing(threshold)
        candidates = values[values .> minimalvalue]
        isempty(candidates) && return missing
        return quantile(candidates, threshold_quantile)
    end

    ismissing(threshold) && return missing
    threshold_value = Float64(threshold)
    isfinite(threshold_value) || return missing
    return threshold_value
end

function _fit_gev_model(xin; min_points::Int=3, fitkwargs...)
    values = _finite_real_values(xin)
    length(values) >= min_points || return missing

    try
        return Extremes.gevfit(values; fitkwargs...)
    catch err
        if err isa MethodError || err isa UndefKeywordError
            rethrow(err)
        end
        return missing
    end
end

function _normalize_extreme_covariate_names(names, keyword::Symbol)
    values = if names isa Symbol || names isa AbstractString
        [Symbol(names)]
    elseif names isa Tuple || names isa AbstractVector
        Symbol.(collect(names))
    else
        throw(ArgumentError("$(keyword) must be a symbol or a collection of symbols."))
    end

    length(unique(values)) == length(values) || throw(ArgumentError("$(keyword) cannot contain duplicate covariate names."))
    return values
end

function _validate_extreme_covariate_values(covariate, covariate_name::Symbol)
    all(value -> ismissing(value) || value isa Real, covariate) || throw(ArgumentError("Covariate $(covariate_name) must contain real values or missing."))
end

function _prepare_extreme_covariates(covariates, selected_names::Vector{Symbol}, ds::YAXArray, fit_dim::Symbol)
    covariates isa NamedTuple || throw(ArgumentError("covariates must be a named tuple."))
    :__climatetools_extreme_value__ in selected_names && throw(ArgumentError("Covariate name `__climatetools_extreme_value__` is reserved."))

    supplied_names = collect(keys(covariates))
    missing_names = setdiff(selected_names, supplied_names)
    isempty(missing_names) || throw(ArgumentError("Missing covariates: $(join(missing_names, ", "))."))

    unused_names = setdiff(supplied_names, selected_names)
    isempty(unused_names) || throw(ArgumentError("Unused covariates: $(join(unused_names, ", "))."))

    response_axis_names = _extreme_axis_names(ds)
    fit_position = findfirst(==(fit_dim), response_axis_names)
    fit_length = size(ds, fit_position)
    prepared = Dict{Symbol, Any}()

    for covariate_name in selected_names
        covariate = getproperty(covariates, covariate_name)

        if covariate isa YAXArray
            covariate_axis_names = _extreme_axis_names(covariate)
            fit_dim in covariate_axis_names || throw(ArgumentError("Covariate $(covariate_name) must contain the fitted dimension $(fit_dim)."))

            extra_axes = setdiff(covariate_axis_names, response_axis_names)
            isempty(extra_axes) || throw(ArgumentError("Covariate $(covariate_name) contains axes not present in the response cube: $(join(extra_axes, ", "))."))

            for axis_name in covariate_axis_names
                collect(lookup(covariate, axis_name)) == collect(lookup(ds, axis_name)) || throw(ArgumentError("Covariate $(covariate_name) coordinates for $(axis_name) must match the response cube."))
            end

            _validate_extreme_covariate_values(covariate, covariate_name)
            prepared[covariate_name] = covariate
        elseif covariate isa AbstractVector
            length(covariate) == fit_length || throw(ArgumentError("Covariate $(covariate_name) length must match the fitted dimension length $(fit_length)."))
            values = collect(covariate)
            _validate_extreme_covariate_values(values, covariate_name)
            prepared[covariate_name] = values
        else
            throw(ArgumentError("Covariate $(covariate_name) must be an AbstractVector or YAXArray."))
        end
    end

    return prepared
end

function _slice_extreme_covariate(covariate::YAXArray, fit_dim::Symbol, keep_axis_names::Tuple, idx::CartesianIndex)
    selectors = Any[Colon() for _ in 1:ndims(covariate)]

    for (position, axis_name) in pairs(_extreme_axis_names(covariate))
        axis_name == fit_dim && continue
        output_position = findfirst(==(axis_name), keep_axis_names)
        selectors[position] = idx[output_position]
    end

    return vec(Array(view(covariate, selectors...)))
end


function _slice_extreme_covariates(prepared::Dict{Symbol, Any}, selected_names::Vector{Symbol}, fit_dim::Symbol, keep_axis_names::Tuple, idx::CartesianIndex)
    sliced = Dict{Symbol, Vector}()

    for covariate_name in selected_names
        covariate = prepared[covariate_name]
        sliced[covariate_name] = covariate isa YAXArray ?
            _slice_extreme_covariate(covariate, fit_dim, keep_axis_names, idx) :
            covariate
    end

    return sliced
end

function _complete_extreme_rows(response, covariates::Dict{Symbol, Vector}, selected_names::Vector{Symbol})
    response_values = collect(response)
    valid = trues(length(response_values))

    for (index, value) in pairs(response_values)
        valid[index] = value isa Real && isfinite(Float64(value))
    end

    for covariate_name in selected_names
        values = covariates[covariate_name]
        for (index, value) in pairs(values)
            valid[index] &= value isa Real && isfinite(Float64(value))
        end
    end

    valid_indices = findall(valid)
    filtered_response = Float64[response_values[index] for index in valid_indices]
    filtered_covariates = Dict(
        covariate_name => Float64[covariates[covariate_name][index] for index in valid_indices]
        for covariate_name in selected_names
    )
    return filtered_response, filtered_covariates, valid_indices
end

function _fit_nonstationary_gev_model(response, covariates::Dict{Symbol, Vector}; locationcovid::Vector{Symbol}, logscalecovid::Vector{Symbol}, shapecovid::Vector{Symbol}, min_points::Int=3, fitkwargs...)
    selected_names = unique(vcat(locationcovid, logscalecovid, shapecovid))
    response_values, covariate_values, valid_indices = _complete_extreme_rows(response, covariates, selected_names)
    parameter_count = 3 + length(locationcovid) + length(logscalecovid) + length(shapecovid)
    length(response_values) >= max(min_points, parameter_count) || return missing, valid_indices

    table = DataFrame(__climatetools_extreme_value__=response_values)
    for covariate_name in selected_names
        table[!, covariate_name] = covariate_values[covariate_name]
    end

    try
        model = Extremes.gevfit(
            table,
            :__climatetools_extreme_value__;
            locationcovid=locationcovid,
            logscalecovid=logscalecovid,
            shapecovid=shapecovid,
            fitkwargs...,
        )
        return model, valid_indices
    catch err
        if err isa MethodError || err isa UndefKeywordError
            rethrow(err)
        end
        return missing, valid_indices
    end
end

function _fit_nonstationary_gev_cube(ds::YAXArray, covariates; dim, locationcovid::Vector{Symbol}, logscalecovid::Vector{Symbol}, shapecovid::Vector{Symbol}, min_points::Int=3, fitkwargs...)
    normalized_dims = _normalize_xmap_dims(dim)
    length(normalized_dims) == 1 || throw(ArgumentError("Non-stationary GEV fitting supports exactly one reduced dimension."))
    _, output_axes, keep_positions, output_size = _extreme_output_layout(ds, normalized_dims)

    fit_dim = only(normalized_dims)
    fit_position = findfirst(==(fit_dim), _extreme_axis_names(ds))
    fit_axis = ds.axes[fit_position]
    selected_names = unique(vcat(locationcovid, logscalecovid, shapecovid))
    prepared_covariates = _prepare_extreme_covariates(covariates, selected_names, ds, fit_dim)
    keep_axis_names = Tuple(_extreme_axis_names(ds)[position] for position in keep_positions)

    models = Array{_EXTREME_MODEL_OUTPUT}(undef, output_size...)
    fill!(models, missing)
    valid_indices = Array{Vector{Int}}(undef, output_size...)

    for idx in CartesianIndices(models)
        model, indices = _fit_nonstationary_gev_model(
            _extreme_slice(ds, keep_positions, idx),
            _slice_extreme_covariates(prepared_covariates, selected_names, fit_dim, keep_axis_names, idx);
            locationcovid=locationcovid,
            logscalecovid=logscalecovid,
            shapecovid=shapecovid,
            min_points=min_points,
            fitkwargs...,
        )
        models[idx] = model
        valid_indices[idx] = indices
    end

    return (
        fit_kind=:gev,
        models=YAXArray(output_axes, models),
        fit_axis=fit_axis,
        valid_indices=YAXArray(output_axes, valid_indices),
    )
end

function _fit_gp_model(xin; threshold=nothing, threshold_quantile::Real=0.95, minimalvalue::Real=1.0, min_exceedances::Int=3, fitkwargs...)
    values = _finite_real_values(xin)
    isempty(values) && return missing, missing

    threshold_value = _gp_threshold(values; threshold=threshold, threshold_quantile=threshold_quantile, minimalvalue=minimalvalue)
    ismissing(threshold_value) && return missing, missing

    exceedances = values[values .> threshold_value] .- threshold_value
    length(exceedances) >= min_exceedances || return missing, threshold_value

    try
        return Extremes.gpfit(exceedances; fitkwargs...), threshold_value
    catch err
        if err isa MethodError || err isa UndefKeywordError
            rethrow(err)
        end
        return missing, threshold_value
    end
end

function _returnlevel_output_layout(model_cube::YAXArray, rlevels)
    level_axis = Dim{:rlevels}(collect(rlevels))
    output_axes = (model_cube.axes..., level_axis)
    output = Array{Union{Missing, Float64}}(undef, size(model_cube)..., length(rlevels))
    fill!(output, missing)
    return output_axes, output
end

function _validate_gev_return_periods(rlevels)
    all(return_period -> return_period isa Real && isfinite(return_period) && return_period > 1, rlevels) || throw(ArgumentError("GEV return periods must be finite real values greater than 1."))
end

function _reject_extreme_covariate_fitkwargs(fitkwargs, fit_kind::Symbol)
    forbidden = fit_kind == :gev ? (:locationcov, :logscalecov, :shapecov) : (:logscalecov, :shapecov)
    supplied = Symbol[keyword for keyword in forbidden if haskey(fitkwargs, keyword)]
    isempty(supplied) || throw(ArgumentError("Raw Extremes.jl covariate keywords are not supported by $(fit_kind)fit_cube: $(join(supplied, ", "))."))
end

function _first_extreme_model(model_cube::YAXArray)
    for model in model_cube
        ismissing(model) || return model
    end

    return nothing
end

_is_gp_model(model) = model isa Extremes.MaximumLikelihoodAbstractExtremeValueModel && model.model isa Extremes.ThresholdExceedance
_is_gev_model(model) = model isa Extremes.MaximumLikelihoodAbstractExtremeValueModel && model.model isa Extremes.BlockMaxima

function _returnlevel_gev_scalar(model, return_period)
    return Extremes.returnlevel(model, return_period).value[1]
end

function _returnlevel_gp_scalar(model, threshold, n_observations::Int, n_obs_per_block::Int, return_period)
    return Extremes.returnlevel(model, threshold, n_observations, n_obs_per_block, return_period).value[1]
end


"""
    gevfit_cube(ds::YAXArray; dim=:time, min_points=3, covariates=nothing,
        locationcovid=Symbol[], logscalecovid=Symbol[], shapecovid=Symbol[],
        fitkwargs...)

Fit a stationary or non-stationary generalized extreme value model at each grid
point by reducing the selected dimension, usually `:time`.

With no selected covariates, the result is a `YAXArray` over the remaining
spatial or ensemble dimensions. Each cell contains either a reusable
`Extremes.MaximumLikelihoodAbstractExtremeValueModel` or `missing` when the
local fit cannot be estimated.

For a non-stationary fit, pass covariate data as a named tuple and select the
covariates used by the location, log-scale, or shape parameter with
`locationcovid`, `logscalecovid`, and `shapecovid`, respectively. These names
and parameter functions follow the `Extremes.gevfit` interface. Each covariate
may be an `AbstractVector` shared by all grid cells or a `YAXArray` containing
the fitted dimension and optionally any subset of the remaining dimensions.
Shared YAXArray axes must have the same coordinates as `ds`.

Non-stationary fitting supports one reduced dimension and returns a named tuple
with fields `fit_kind`, `models`, `fit_axis`, and `valid_indices`. Missing and
non-finite response or covariate rows are removed jointly at each grid cell;
`valid_indices` records their positions on `fit_axis`. A fit requires at least
`min_points` complete rows and at least as many complete rows as regression
coefficients. Model field `θ̂` contains the fitted regression coefficients,
while `Extremes.params(model)` evaluates the GEV location, scale, and shape at
each retained covariate row. Pass the fit bundle to `returnlevel_cube` to obtain
effective return levels on axes `(remaining axes..., fit_axis, :rlevels)`.

`gevfit_cube` operates directly on block maxima and therefore does not accept
threshold-filtering keywords such as `minimalvalue`.
"""
function gevfit_cube(ds::YAXArray; dim=:time, min_points::Int=3, covariates=nothing, locationcovid=Symbol[], logscalecovid=Symbol[], shapecovid=Symbol[], fitkwargs...)
    if haskey(fitkwargs, :minimalvalue)
        throw(ArgumentError("gevfit_cube does not accept `minimalvalue`; pass block maxima directly and use gpfit_cube for threshold exceedance models."))
    end
    _reject_extreme_covariate_fitkwargs(fitkwargs, :gev)

    location_names = _normalize_extreme_covariate_names(locationcovid, :locationcovid)
    logscale_names = _normalize_extreme_covariate_names(logscalecovid, :logscalecovid)
    shape_names = _normalize_extreme_covariate_names(shapecovid, :shapecovid)
    selected_names = unique(vcat(location_names, logscale_names, shape_names))

    if !isempty(selected_names)
        isnothing(covariates) && throw(ArgumentError("covariates must be provided for a non-stationary GEV fit."))
        return _fit_nonstationary_gev_cube(
            ds,
            covariates;
            dim=dim,
            locationcovid=location_names,
            logscalecovid=logscale_names,
            shapecovid=shape_names,
            min_points=min_points,
            fitkwargs...,
        )
    end

    if !isnothing(covariates)
        covariates isa NamedTuple || throw(ArgumentError("covariates must be a named tuple."))
        isempty(keys(covariates)) || throw(ArgumentError("At least one of locationcovid, logscalecovid, or shapecovid must select the supplied covariates."))
    end

    return _fit_extreme_model_cube(ds; reduced_dims=dim) do xin
        _fit_gev_model(xin; min_points=min_points, fitkwargs...)
    end
end


"""
    gpfit_cube(ds::YAXArray; dim=:time, threshold=nothing, threshold_quantile=0.95,
        minimalvalue=1.0, min_exceedances=3, return_thresholds=false,
        n_obs_per_block=365, fitkwargs...)

Fit a stationary generalized Pareto model at each grid point by reducing the
selected dimension, usually `:time`.

When `threshold=nothing`, a local threshold is estimated independently at each
grid point from the `threshold_quantile` of the values above `minimalvalue`
(for example `0.95`, i.e. the 95th percentile).
The fitted model is estimated on threshold exceedances, so `return_thresholds=true`
returns a named tuple that can be reused later for return-level calculations on
the original scale.
"""
function gpfit_cube(ds::YAXArray; dim=:time, threshold=nothing, threshold_quantile::Real=0.95, minimalvalue::Real=1.0, min_exceedances::Int=3, return_thresholds::Bool=false, n_obs_per_block::Int=365, fitkwargs...)
    _reject_extreme_covariate_fitkwargs(fitkwargs, :gp)
    _, output_axes, keep_positions, output_size = _extreme_output_layout(ds, dim)
    models = Array{_EXTREME_MODEL_OUTPUT}(undef, output_size...)
    fill!(models, missing)

    thresholds = nothing
    if return_thresholds
        thresholds = Array{Union{Missing, Float64}}(undef, output_size...)
        fill!(thresholds, missing)
    end

    for idx in CartesianIndices(models)
        model, threshold_value = _fit_gp_model(
            _extreme_slice(ds, keep_positions, idx);
            threshold=threshold,
            threshold_quantile=threshold_quantile,
            minimalvalue=minimalvalue,
            min_exceedances=min_exceedances,
            fitkwargs...,
        )
        models[idx] = model

        if return_thresholds
            thresholds[idx] = threshold_value
        end
    end

    model_cube = YAXArray(output_axes, models)
    return return_thresholds ? (
        models=model_cube,
        thresholds=YAXArray(output_axes, thresholds),
        n_observations=_extreme_n_observations(ds, dim),
        n_obs_per_block=n_obs_per_block,
    ) : model_cube
end

"""
    returnlevel_cube(model_cube::YAXArray; rlevels=[2, 5, 10, 25, 50, 100, 1000])

Compute return levels from a cube of stationary GEV models returned by
`gevfit_cube`. The output appends a `:rlevels` axis after the remaining
dimensions.
"""
function returnlevel_cube(model_cube::YAXArray; rlevels=[2, 5, 10, 25, 50, 100, 1000])
    _validate_gev_return_periods(rlevels)
    first_model = _first_extreme_model(model_cube)
    output_axes, output = _returnlevel_output_layout(model_cube, rlevels)

    isnothing(first_model) && return YAXArray(output_axes, output)

    if _is_gp_model(first_model)
        error("GP model cubes require thresholds and observation metadata. Pass the named tuple returned by gpfit_cube(...; return_thresholds=true) or call returnlevel_cube(models, thresholds; n_observations=..., n_obs_per_block=...).")
    end

    _is_gev_model(first_model) || error("Unsupported model type in returnlevel_cube.")

    model_array = Array(model_cube)
    for idx in CartesianIndices(model_array)
        model = model_array[idx]
        ismissing(model) && continue

        for (level_index, return_period) in pairs(rlevels)
            output[(Tuple(idx)..., level_index)...] = _returnlevel_gev_scalar(model, return_period)
        end
    end

    return YAXArray(output_axes, output)
end

"""
    returnlevel_cube(models::YAXArray, thresholds::YAXArray;
        n_observations, n_obs_per_block=365,
        rlevels=[1, 2, 5, 10, 25, 50, 100, 1000])

Compute return levels from a cube of stationary GP models and the companion
threshold cube used to fit them.
"""
function returnlevel_cube(models::YAXArray, thresholds::YAXArray; n_observations::Int, n_obs_per_block::Int=365, rlevels=[1, 2, 5, 10, 25, 50, 100, 1000])
    size(models) == size(thresholds) || error("Model and threshold cubes must have the same size.")

    output_axes, output = _returnlevel_output_layout(models, rlevels)
    model_array = Array(models)
    threshold_array = Array(thresholds)

    for idx in CartesianIndices(model_array)
        model = model_array[idx]
        threshold = threshold_array[idx]
        (ismissing(model) || ismissing(threshold)) && continue

        for (level_index, return_period) in pairs(rlevels)
            output[(Tuple(idx)..., level_index)...] = _returnlevel_gp_scalar(
                model,
                threshold,
                n_observations,
                n_obs_per_block,
                return_period,
            )
        end
    end

    return YAXArray(output_axes, output)
end

function _returnlevel_nonstationary_gev_cube(fit_bundle::NamedTuple; rlevels)
    _validate_gev_return_periods(rlevels)
    hasproperty(fit_bundle, :models) || error("Non-stationary GEV fit bundle must contain a `models` field.")
    hasproperty(fit_bundle, :fit_axis) || error("Non-stationary GEV fit bundle must contain a `fit_axis` field.")
    hasproperty(fit_bundle, :valid_indices) || error("Non-stationary GEV fit bundle must contain a `valid_indices` field.")

    models = fit_bundle.models
    valid_indices = fit_bundle.valid_indices
    models isa YAXArray || error("Non-stationary GEV fit bundle `models` must be a YAXArray.")
    valid_indices isa YAXArray || error("Non-stationary GEV fit bundle `valid_indices` must be a YAXArray.")
    size(models) == size(valid_indices) || error("Non-stationary GEV model and valid-index cubes must have the same size.")
    _extreme_axis_names(models) == _extreme_axis_names(valid_indices) || error("Non-stationary GEV model and valid-index cubes must have the same axes.")
    for axis_name in _extreme_axis_names(models)
        collect(lookup(models, axis_name)) == collect(lookup(valid_indices, axis_name)) || error("Non-stationary GEV model and valid-index cubes must have the same axis coordinates.")
    end

    fit_axis = fit_bundle.fit_axis
    output_axes = (models.axes..., fit_axis, Dim{:rlevels}(collect(rlevels)))
    output = Array{Union{Missing, Float64}}(undef, size(models)..., length(fit_axis), length(rlevels))
    fill!(output, missing)

    model_array = Array(models)
    valid_index_array = Array(valid_indices)

    for idx in CartesianIndices(model_array)
        model = model_array[idx]
        ismissing(model) && continue
        _is_gev_model(model) || error("Non-stationary GEV fit bundle contains an unsupported model type.")

        rows = valid_index_array[idx]
        all(row -> 1 <= row <= length(fit_axis), rows) || error("Non-stationary GEV valid indices must reference the fitted axis.")
        length(unique(rows)) == length(rows) || error("Non-stationary GEV valid indices must be unique.")

        for (level_index, return_period) in pairs(rlevels)
            values = Extremes.returnlevel(model, return_period).value
            length(values) == length(rows) || error("Non-stationary GEV return-level values do not match the stored valid indices.")

            for (value_index, row) in pairs(rows)
                output[(Tuple(idx)..., row, level_index)...] = values[value_index]
            end
        end
    end

    return YAXArray(output_axes, output)
end

"""
    returnlevel_cube(fit_bundle::NamedTuple; rlevels=nothing)

Compute return levels from a non-stationary GEV fit bundle returned by
`gevfit_cube(...; covariates=...)` or a GP fit bundle returned by
`gpfit_cube(...; return_thresholds=true)`.
"""
function returnlevel_cube(fit_bundle::NamedTuple; rlevels=nothing)
    if hasproperty(fit_bundle, :fit_kind) && fit_bundle.fit_kind == :gev
        gev_rlevels = isnothing(rlevels) ? [2, 5, 10, 25, 50, 100, 1000] : rlevels
        return _returnlevel_nonstationary_gev_cube(fit_bundle; rlevels=gev_rlevels)
    end

    hasproperty(fit_bundle, :models) || error("GP fit bundle must contain a `models` field.")
    hasproperty(fit_bundle, :thresholds) || error("GP fit bundle must contain a `thresholds` field.")
    hasproperty(fit_bundle, :n_observations) || error("GP fit bundle must contain an `n_observations` field.")
    hasproperty(fit_bundle, :n_obs_per_block) || error("GP fit bundle must contain an `n_obs_per_block` field.")
    gp_rlevels = isnothing(rlevels) ? [1, 2, 5, 10, 25, 50, 100, 1000] : rlevels

    return returnlevel_cube(
        fit_bundle.models,
        fit_bundle.thresholds;
        n_observations=fit_bundle.n_observations,
        n_obs_per_block=fit_bundle.n_obs_per_block,
        rlevels=gp_rlevels,
    )
end
