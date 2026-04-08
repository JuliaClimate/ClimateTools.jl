function rlevels_cube(xout, xin; threshold=nothing, rlevels = [1, 2, 5, 10, 25, 50, 100, 1000], minimalvalue=1.0)

    if all(ismissing, xin)
        xout .= missing
        return    
    end

    if isnothing(threshold)
        dataforquantile = xin[xin .> minimalvalue]
        if  !isempty(dataforquantile)
            threshold = quantile(dataforquantile, 0.92)
            exceedances = xin[xin .> threshold] .- threshold
            model = Extremes.gpfit(exceedances)

            nobs = size(xin,1)        

            nobsperblock = 365

            r_h = returnlevel.(model, threshold, nobs, nobsperblock, rlevels)

            xout .= [x.value[1] for x in r_h]

        else       
            xout .= missing
            return
        end
        
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
    catch
        return missing
    end
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
    catch
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
    gevfit_cube(ds::YAXArray; dim=:time, min_points=3, fitkwargs...)

Fit a stationary generalized extreme value model at each grid point by reducing
the selected dimension, usually `:time`.

The result is a `YAXArray` over the remaining spatial or ensemble dimensions,
with each cell containing either a reusable `Extremes.MaximumLikelihoodAbstractExtremeValueModel`
or `missing` when the local fit cannot be estimated.
"""
function gevfit_cube(ds::YAXArray; dim=:time, min_points::Int=3, fitkwargs...)
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
    returnlevel_cube(model_cube::YAXArray; rlevels=[1, 2, 5, 10, 25, 50, 100, 1000])

Compute return levels from a cube of stationary GEV models returned by
`gevfit_cube`. The output appends a `:rlevels` axis after the remaining
dimensions.
"""
function returnlevel_cube(model_cube::YAXArray; rlevels=[1, 2, 5, 10, 25, 50, 100, 1000])
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

"""
    returnlevel_cube(fit_bundle::NamedTuple; rlevels=[1, 2, 5, 10, 25, 50, 100, 1000])

Compute return levels from the bundle returned by
`gpfit_cube(...; return_thresholds=true)`.
"""
function returnlevel_cube(fit_bundle::NamedTuple; rlevels=[1, 2, 5, 10, 25, 50, 100, 1000])
    hasproperty(fit_bundle, :models) || error("GP fit bundle must contain a `models` field.")
    hasproperty(fit_bundle, :thresholds) || error("GP fit bundle must contain a `thresholds` field.")
    hasproperty(fit_bundle, :n_observations) || error("GP fit bundle must contain an `n_observations` field.")
    hasproperty(fit_bundle, :n_obs_per_block) || error("GP fit bundle must contain an `n_obs_per_block` field.")

    return returnlevel_cube(
        fit_bundle.models,
        fit_bundle.thresholds;
        n_observations=fit_bundle.n_observations,
        n_obs_per_block=fit_bundle.n_obs_per_block,
        rlevels=rlevels,
    )
end
