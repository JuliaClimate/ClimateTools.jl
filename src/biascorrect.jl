function _time_dim_symbol(ds)
    try
        lookup(ds, :time)
        return :time
    catch
    end
    try
        lookup(ds, :Ti)
        return :Ti
    catch
    end
    error("Could not infer time dimension symbol. Expected :time.")
end

"""
    qqmap(obs::YAXArray, ref::YAXArray, fut::YAXArray; method::String="Additive", detrend::Bool=true, order::Int=4, window::Int=15, rankn::Int=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, wet_threshold::Real=1.0, interp=Linear(), extrap=Interpolations.Flat())

This function performs quantile mapping bias correction on data.

## Arguments
- `obs::YAXArray`: The observed data.
- `ref::YAXArray`: The reference data.
- `fut::YAXArray`: The future data.
- `method::String`: The method to apply on raw data for bias correction. Default is "Additive".
- `detrend::Bool`: Whether to detrend the data before bias correction. Default is `true`.
- `order::Int`: The order of the polynomial used for detrending. Default is 4.
- `window::Int`: The window size (+/- window relative to a given julian day) for calculating the Julian days. Default is 15.
- `rankn::Int`: The number of quantiles to use for mapping. Default is 50.
- `qmin::Real`: The minimum quantile value. Default is 0.01.
- `qmax::Real`: The maximum quantile value. Default is 0.99.
- `thresnan::Float64`: The threshold for bias correcting a grid based in presence of NaN values. Default is 0.1 (i.e. if there is more than 10% of NaN values, the grid point is not corrected).
- `keep_original::Bool`: Whether to keep the original data in the output if there is more than the threshold of NaN values. Default is `false`.
- `wet_threshold::Real`: Minimum precipitation amount used to build the multiplicative transfer for nonnegative series. Values below this threshold keep their original intensity. Default is `1.0`.
- `interp`: The interpolation method to use between quantiles. Default is `Linear()`.
- `extrap`: The extrapolation method to use over qmin and qmax. Default is `Interpolations.Flat()`.

## Returns
- A bias-corrected YAXArray.

# References

- Themeßl, M. J., Gobiet, A., and Leuprecht, A. (2012). Empirical-statistical downscaling and error correction of regional climate models and its impact on the climate change signal.
- Piani, C., Haerter, J. O., and Coppola, E. (2010). Statistical bias correction for daily precipitation in regional climate models over Europe.
- Roy, P., Rondeau-Genesse, G., Jalbert, J., and Fournier, E. (2023). Climate scenarios of extreme precipitation using a combination of parametric and non-parametric bias correction methods in the province of Québec. DOI: 10.1080/07011784.2023.2220682.

"""
function qqmap(obs::YAXArray, ref::YAXArray, fut::YAXArray; method::String="Additive", detrend::Bool=true, order::Int=4, window::Int=15, rankn::Int=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, wet_threshold::Real=1.0, interp=Linear(), extrap=Interpolations.Flat())
    timedim = _time_dim_symbol(obs)
    
    # Aligning calendars
    obs = drop29thfeb(obs, dimt=timedim)
    ref = drop29thfeb(ref, dimt=timedim)
    fut = drop29thfeb(fut, dimt=timedim)

    # Julian days vectors
    obs_jul = Dates.dayofyear.(lookup(obs, timedim))
    ref_jul = Dates.dayofyear.(lookup(ref, timedim))
    fut_jul = Dates.dayofyear.(lookup(fut, timedim))
    
    # # Monthly values of time vector
    # obs_month = Dates.month.(lookup(obs, :time))
    # ref_month = Dates.month.(lookup(ref, :time))
    # fut_month = Dates.month.(lookup(fut, :time))
    
    # # Seasonal values of time vector
    # obs_season = get_seasons(obs_month)
    # ref_season = get_seasons(ref_month)
    # fut_season = get_seasons(fut_month)
    
    
    if Base.minimum(ref_jul) == 1 && Base.maximum(ref_jul) == 365
        days = 1:365
    else
        days = Base.minimum(ref_jul) + window:Base.maximum(ref_jul) - window
    end
    
    indims = InDims(string(timedim))
    outdims = OutDims(Dim{:time}(lookup(fut, timedim)))

    return mapCube(qqmap, (obs, ref, fut), indims=(indims, indims, indims), outdims=outdims, days, obs_jul, ref_jul, fut_jul, method=method, detrend=detrend, order=order, window=window, rankn=rankn, qmin=qmin, qmax=qmax, thresnan=thresnan, keep_original=keep_original, wet_threshold=wet_threshold, interp=interp, extrap=extrap, nthreads=Threads.nthreads())
    
end

"""
    qqmap(dataout, obsvec, refvec, futvec, days, obs_jul, ref_jul, fut_jul; method::String="Additive", detrend::Bool=true, order::Int=4, window::Int64=15, rankn::Int64=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, wet_threshold::Real=1.0, interp=Linear(), extrap=Interpolations.Flat())

In-place kernel used by `qqmap` for a single grid-point time series.

It corrects the values in `futvec` for each day specified in `days` based on the quantiles
estimated from the corresponding moving-window samples in `obsvec` and `refvec`.

## Arguments
- `dataout`: Output array where the corrected values will be stored.
- `obsvec`: Array of observed values.
- `refvec`: Array of reference values.
- `futvec`: Array of future values to be corrected.
- `days`: Array of days of the year for which correction will be performed.
- `obs_jul`: Array of Julian day values corresponding to `obsvec`.
- `ref_jul`: Array of Julian day values corresponding to `refvec`.
- `fut_jul`: Array of Julian day values corresponding to `futvec`.

## Optional Arguments
- `method::String`: Method used for correction. Default is "Additive".
- `detrend::Bool`: Whether to detrend the data before correction. Default is `true`.
- `order::Int`: Order of the polynomial used for detrending. Default is 4.
- `window::Int64`: Size of the moving window for selecting reference values. Default is 15.
- `rankn::Int64`: Number of quantiles to estimate. Default is 50.
- `qmin::Real`: Minimum quantile value. Default is 0.01.
- `qmax::Real`: Maximum quantile value. Default is 0.99.
- `thresnan::Float64`: Threshold for the percentage of NaN values allowed in the data. Default is 0.1.
- `keep_original::Bool`: Whether to keep the original values in case of too many NaN values. Default is `false`.
- `wet_threshold::Real`: Minimum precipitation amount used to build multiplicative transfer functions for nonnegative series. Default is `1.0`.
- `interp`: Interpolation method used for building the correction function. Default is `Linear()`.
- `extrap`: Extrapolation method used for building the correction function. Default is `Interpolations.Flat()`.

The function modifies `dataout` in-place and returns nothing.
"""
function qqmap(dataout, obsvec, refvec, futvec, days, obs_jul, ref_jul, fut_jul; method::String="Additive", detrend::Bool=true, order::Int=4, window::Int64=15, rankn::Int64=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, wet_threshold::Real=1.0, interp=Linear(), extrap=Interpolations.Flat())
    

    # range over which quantiles are estimated
    P = range(qmin, stop=qmax, length = rankn)
    obsP = similar(obsvec, length(P))
    refP = similar(refvec, length(P))
    sf_refP = similar(refvec, length(P))
    
    if detrend == true
        obsvec, _ = _detrend_series(obsvec, order=order)
        refvec, _ = _detrend_series(refvec, order=order)
        futvec, poly_values = _detrend_series(futvec, order=order)
    end

    # LOOP OVER ALL DAYS OF THE YEAR
    for ijulian in days

        # idx for values we want to correct
        idxfut = (fut_jul .== ijulian)

        # Find all index of moving window around ijulian day of year
        idxobs = find_julianday_idx(obs_jul, ijulian, window)
        idxref = find_julianday_idx(ref_jul, ijulian, window)

        obsval = obsvec[idxobs]
        refval = refvec[idxref]
        futval = futvec[idxfut]
        
        obsval_unique = unique(obsval)
        refval_unique = unique(refval)
        futval_unique = unique(futval)

        # Random perturbation is added to avoid division by zero and duplicate values
        if lowercase(method) == "multiplicative"
            # Random vectors parameters
            init     = 0.000001
            interval = 0.00001
            stop     = 0.001      

            # generate_random_perturbations!(obsval, obsP)#init, interval, stop)
            # generate_random_perturbations!(refval, refP)#, interval, stop)
            # generate_random_perturbations!(futval, init, interval, stop)
            
        end        

        obsvalnan = isnan.(obsval)
        refvalnan = isnan.(refval)
        futvalnan = isnan.(futval)

        if (sum(obsvalnan) < (length(obsval) * thresnan)) & (sum(refvalnan) < (length(refval) * thresnan)) & (sum(futvalnan) < (length(futval) * thresnan))

            # Estimate quantiles for obs and ref for ijulian
            quantile!(obsP, obsval[.!obsvalnan], P, sorted = false)
            quantile!(refP, refval[.!refvalnan], P, sorted = false)

            if lowercase(method) == "additive" # used for temperature
                sf_refP .= obsP .- refP
                
                # Build interpolation
                itp = build_itp(refP, sf_refP, interp, extrap)
                
                # dataout[idxfut] .= itp(futval) .+ futval
                futval .= itp(futval) .+ futval

            elseif lowercase(method) == "multiplicative" # used for precipitation
                _apply_multiplicative_qqmap!(futval, obsP, refP, sf_refP, interp, extrap;
                    wet_threshold=wet_threshold,
                    wet_day_guard=_should_use_wet_day_guard(obsval[.!obsvalnan], refval[.!refvalnan], futval[.!futvalnan]),
                )

            else
                error("Wrong method")
            end
            # Replace values with new ones
            dataout[idxfut] .= futval
        else

            if keep_original
                # # Replace values with original ones (i.e. too may NaN values for robust quantile estimation)
                dataout[idxfut] .= futval
            else
                # DO NOTHING (e.g. if there is no reference, we want NaNs and not original values)
                dataout[idxfut] .= NaN
            end
        end
    end
    
    if detrend == true
        dataout .= dataout .+ poly_values
    end

end

"""
    qqmap_bulk(obs::YAXArray, ref::YAXArray, fut::YAXArray; method::String="Additive", detrend::Bool=true, order::Int=4, window::Int=15, rankn::Int=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, wet_threshold::Real=1.0, interp=Linear(), extrap=Interpolations.Flat())

This function performs bulk quantile mapping bias correction on data.

Unlike `qqmap`, this method estimates a single transfer function over the full sample rather than
building day-of-year-specific mappings with a moving seasonal window.

## Arguments
- `obs::YAXArray`: The observed data.
- `ref::YAXArray`: The reference data.
- `fut::YAXArray`: The future data.
- `method::String`: The method to apply on raw data for bias correction. Default is "Additive".
- `detrend::Bool`: Whether to detrend the data before bias correction. Default is `true`.
- `order::Int`: The order of the polynomial used for detrending. Default is 4.
- `rankn::Int`: The number of quantiles to use for mapping. Default is 50.
- `qmin::Real`: The minimum quantile value. Default is 0.01.
- `qmax::Real`: The maximum quantile value. Default is 0.99.
- `thresnan::Float64`: The threshold for bias correcting a grid based in presence of NaN values. Default is 0.1.
- `keep_original::Bool`: Whether to keep the original data in the output if there is more than the threshold of NaN values. Default is `false`.
- `wet_threshold::Real`: Minimum precipitation amount used to build the multiplicative transfer for nonnegative series. Values below this threshold keep their original intensity. Default is `1.0`.
- `interp`: The interpolation method to use between quantiles. Default is `Linear()`.
- `extrap`: The extrapolation method to use over qmin and qmax. Default is `Interpolations.Flat()`.

## Returns
- A bias-corrected YAXArray.

# References

- Themeßl, M. J., Gobiet, A., and Leuprecht, A. (2012). Empirical-statistical downscaling and error correction of regional climate models and its impact on the climate change signal.
- Piani, C., Haerter, J. O., and Coppola, E. (2010). Statistical bias correction for daily precipitation in regional climate models over Europe.
- Roy, P., Rondeau-Genesse, G., Jalbert, J., and Fournier, E. (2023). Climate scenarios of extreme precipitation using a combination of parametric and non-parametric bias correction methods in the province of Québec. DOI: 10.1080/07011784.2023.2220682.
"""
function qqmap_bulk(obs::YAXArray, ref::YAXArray, fut::YAXArray; method::String="Additive", detrend::Bool=true, order::Int=4, window::Int=15, rankn::Int=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, wet_threshold::Real=1.0, interp=Linear(), extrap=Interpolations.Flat())
    timedim = _time_dim_symbol(obs)


    # Aligning calendars
    obs = drop29thfeb(obs, dimt=timedim)
    ref = drop29thfeb(ref, dimt=timedim)
    fut = drop29thfeb(fut, dimt=timedim)
    
    indims = InDims(string(timedim))
    outdims = OutDims(Dim{:time}(lookup(fut, timedim)))

    return mapCube(qqmap_bulk, (obs, ref, fut), indims=(indims, indims, indims), outdims=outdims, method=method, detrend=detrend, order=order, rankn=rankn, qmin=qmin, qmax=qmax, thresnan=thresnan, keep_original=keep_original, wet_threshold=wet_threshold, interp=interp, extrap=extrap, nthreads=Threads.nthreads())
    
end


function qqmap_bulk(dataout, obsvec, refvec, futvec; method::String="Additive", detrend::Bool=true, order::Int=4, rankn::Int64=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, wet_threshold::Real=1.0, interp=Linear(), extrap=Interpolations.Flat())
    

    # range over which quantiles are estimated
    P = range(qmin, stop=qmax, length = rankn)
    obsP = similar(obsvec, length(P))
    refP = similar(refvec, length(P))
    sf_refP = similar(refvec, length(P))

    if detrend == true
        obsvec, _ = _detrend_series(obsvec, order=order)
        refvec, _ = _detrend_series(refvec, order=order)
        futvec, poly_values = _detrend_series(futvec, order=order)
    end

    obsvecnan = isnan.(obsvec)
    refvecnan = isnan.(refvec)
    futvecnan = isnan.(futvec)

    if (sum(obsvecnan) < (length(obsvec) * thresnan)) & (sum(refvecnan) < (length(refvec) * thresnan)) & (sum(futvecnan) < (length(futvec) * thresnan))

        # Estimate quantiles for obs and ref for ijulian
        quantile!(obsP, obsvec[.!obsvecnan], P, sorted = false)
        quantile!(refP, refvec[.!refvecnan], P, sorted = false)

        if lowercase(method) == "additive" # used for temperature
            sf_refP .= obsP .- refP
            
            # Build interpolation
            itp = build_itp(refP, sf_refP, interp, extrap)
            
            # dataout[idxfut] .= itp(futval) .+ futval
            futvec .= itp(futvec) .+ futvec

        elseif lowercase(method) == "multiplicative" # used for precipitation
            _apply_multiplicative_qqmap!(futvec, obsP, refP, sf_refP, interp, extrap;
                wet_threshold=wet_threshold,
                wet_day_guard=_should_use_wet_day_guard(obsvec[.!obsvecnan], refvec[.!refvecnan], futvec[.!futvecnan]),
            )

        else
            error("Wrong method")
        end
        # Replace values with new ones
        dataout .= futvec
    else

        if keep_original
            # # Replace values with original ones (i.e. too may NaN values for robust quantile estimation)
            dataout .= futvec
        else
            # DO NOTHING (e.g. if there is no reference, we want NaNs and not original values)
            dataout .= NaN
        end
    end

    if detrend == true
        dataout .= dataout .+ poly_values
    end

end

function _is_nonnegative_series(vec; atol::Real=1e-8)
    finite = vec[isfinite.(vec)]
    isempty(finite) && return false
    return Base.minimum(finite) >= -atol
end

function _should_use_wet_day_guard(obsvec, refvec, futvec; wet_threshold::Real=1.0)
    wet_threshold <= 0 && return false
    return _is_nonnegative_series(obsvec) && _is_nonnegative_series(refvec) && _is_nonnegative_series(futvec)
end

function _apply_multiplicative_qqmap!(futvec, obsP, refP, sf_refP, interp, extrap; wet_threshold::Real=1.0, wet_day_guard::Bool=false)
    sf_refP .= obsP ./ refP
    sf_refP[sf_refP .< 0] .= 0.0
    sf_refP[isnan.(sf_refP)] .= 0.0
    sf_refP[isinf.(sf_refP)] .= 0.0

    if wet_day_guard
        wetmask = refP .>= wet_threshold

        if count(wetmask) < 2
            futvec .= max.(futvec, 0.0)
            return futvec
        end

        ref_knots = copy(refP[wetmask])
        scale_knots = copy(sf_refP[wetmask])
        itp = build_itp(ref_knots, scale_knots, interp, extrap)

        wetfut = futvec .>= wet_threshold
        futvec[wetfut] .= itp(futvec[wetfut]) .* futvec[wetfut]
        futvec[.!wetfut] .= max.(futvec[.!wetfut], 0.0)
        return futvec
    end

    itp = build_itp(copy(refP), copy(sf_refP), interp, extrap)
    futvec .= itp(futvec) .* futvec
    return futvec
end

function _find_dim_in_axes(cube::YAXArray, candidates::Tuple)
    names = name.(cube.axes)
    for cand in candidates
        if cand in names
            return cand
        end
    end
    return nothing
end

function _flatten_coord_values(cube::YAXArray, spatial_dims::Vector{Symbol}, target_dim::Symbol)
    target_pos = findfirst(==(target_dim), spatial_dims)
    target_pos === nothing && error("Could not find dimension $(target_dim) in spatial dimensions.")

    spatial_sizes = Tuple(length(lookup(cube, d)) for d in spatial_dims)
    target_lookup = lookup(cube, target_dim)
    out = Vector{Float64}(undef, prod(spatial_sizes))

    for (i, idx) in enumerate(CartesianIndices(spatial_sizes))
        out[i] = Float64(target_lookup[idx[target_pos]])
    end

    return out
end

function _cube_to_timeseries_matrix(cube::YAXArray, timedim::Symbol, spatial_dims::Vector{Symbol})
    names = collect(name.(cube.axes))
    perm = Int[]

    for d in spatial_dims
        idx = findfirst(==(d), names)
        idx === nothing && error("Missing dimension $(d) in input cube.")
        push!(perm, idx)
    end

    tidx = findfirst(==(timedim), names)
    tidx === nothing && error("Missing time dimension $(timedim) in input cube.")
    push!(perm, tidx)

    arr = permutedims(Array(cube), Tuple(perm))
    spatial_shape = size(arr)[1:end-1]
    mat = reshape(arr, :, size(arr, ndims(arr)))

    return mat, spatial_shape
end

function _matrix_to_cube_array(data::AbstractMatrix, spatial_shape, spatial_dims::Vector{Symbol}, timedim::Symbol, target_names::Vector{Symbol})
    reshaped = reshape(data, (spatial_shape..., size(data, 2)))
    current_names = vcat(spatial_dims, [timedim])
    perm = [findfirst(==(n), current_names) for n in target_names]
    any(isnothing.(perm)) && error("Could not map current dimensions back to target dimensions.")

    return permutedims(reshaped, Tuple(Int.(perm)))
end

function _get_max_clusters(clusters)
    out = Float64[]
    for c in clusters
        append!(out, c.value)
    end
    return out
end

function _get_position_clusters(clusters)
    out = Int[]
    for c in clusters
        append!(out, c.position)
    end
    return out
end

function _fit_gpd_distribution(vec::Vector{Float64}, thres::Real; runlength::Int=2)
    try
        clusters = Extremes.getcluster(vec, thres, runlength=runlength)
        maxvalues = _get_max_clusters(clusters)
        positions = _get_position_clusters(clusters)
        isempty(maxvalues) && return nothing, Int[], Float64[]

        exceedances = maxvalues .- thres
        model = Extremes.gpfit(exceedances)
        dist = Extremes.getdistribution(model)[1]
        return dist, positions, maxvalues
    catch
        return nothing, Int[], Float64[]
    end
end

function _threshold_from_vectors(obsvec::Vector{Float64}, refvec::Vector{Float64}; P::Real=0.95)
    obsvals = [x for x in obsvec if isfinite(x) && x >= 1.0]
    refvals = [x for x in refvec if isfinite(x) && x >= 1.0]
    if isempty(obsvals) || isempty(refvals)
        return NaN
    end
    return mean((quantile(obsvals, P), quantile(refvals, P)))
end

function gev2gpd(mu::Real, sigma::Real, xi::Real, thres::Real)
    return sigma + xi * (thres - mu)
end

function _extremes_weight(x; frac=0.25, power=1.0)
    denom = (Base.maximum(x) - Base.minimum(x)) * frac
    if denom == 0
        return ones(length(x))
    end
    w = ((x .- Base.minimum(x)) ./ denom) .^ power
    w[w .> 1.0] .= 1.0
    return w
end

function estimate_gpd(lat::Real, lon::Real, gevparams::DataFrame, thres::Real)
    lats = Float64.(Array(gevparams[!, :lat]))
    lons = Float64.(Array(gevparams[!, :lon]))
    d2 = (lats .- Float64(lat)) .^ 2 .+ (lons .- Float64(lon)) .^ 2
    idx = Base.argmin(d2)

    mu = gevparams[idx, :mu]
    sigma = gevparams[idx, :sigma]
    xi = gevparams[idx, :xi]
    sigma2 = gev2gpd(mu, sigma, xi, thres)

    return Extremes.GeneralizedPareto(sigma2, sigma, xi)
end

function build_idx_biascorrect(refvec::Vector{Float64}, futvec::Vector{Float64})
    w_length = length(refvec)
    fut_l = length(futvec)

    R = w_length / 2:w_length / 3:fut_l

    idxi = Vector{Int}(undef, length(R))
    idxf = Vector{Int}(undef, length(R))
    idxi_corr = Vector{Int}(undef, length(R))
    idxf_corr = Vector{Int}(undef, length(R))

    for c in eachindex(R)
        idxi[c] = round(Int, max(R[c] - w_length / 2 + 1, 1))
        idxf[c] = round(Int, min(R[c] + w_length / 2, fut_l))

        idxi_corr[c] = round(Int, idxi[c] + (w_length / 2 - (w_length / 2) / 3) - 1)
        if c == 1
            idxi_corr[c] = 1
        end

        if c == length(R)
            idxf_corr[c] = max(round(Int, idxi[c] + w_length - w_length / 3), fut_l)
        else
            idxf_corr[c] = round(Int, idxi[c] + (w_length / 2 + (w_length / 2) / 3) - 1)
        end
    end

    if length(idxi_corr) > 1
        idxi_corr[2:end] .-= 1
    end

    return idxi, idxf, idxi_corr, idxf_corr
end

"""
    biascorrect_extremes(obs::YAXArray, ref::YAXArray, fut::YAXArray; kwargs...)

Bias-correct future data with multiplicative QQ mapping and an extreme-tail GPD adjustment.
The base correction is done by `qqmap(..., method="multiplicative")`, then identified
cluster extremes in `fut` are remapped from the fitted future tail distribution to either
an observed-tail GPD (estimated from `obs`) or an externally provided parameter table.

This function implements the mixed QQM-GPD strategy described by Roy et al. (2023),
"Climate scenarios of extreme precipitation using a combination of parametric and
non-parametric bias correction methods in the province of Québec"
(DOI: 10.1080/07011784.2023.2220682).

In the paper's notation, the bulk term `H_X` corresponds to the multiplicative `qqmap`
baseline used here, the tail term `G` corresponds to the GPD models produced by
`_fit_gpd_distribution` or `estimate_gpd`, the lower threshold `ℓ` corresponds to the
threshold computed by `_threshold_from_vectors` and controlled by `P`, the upper
transition threshold `t` is represented operationally through `_extremes_weight` via
`frac` and `power`, and the transition weight `w` is the blend returned by
`_extremes_weight` when tail-corrected values are merged back into the `qqmap` base
field.

# References

- Roy, P., Rondeau-Genesse, G., Jalbert, J., and Fournier, E. (2023). Climate scenarios of extreme precipitation using a combination of parametric and non-parametric bias correction methods in the province of Québec. DOI: 10.1080/07011784.2023.2220682.
"""
function biascorrect_extremes(obs::YAXArray, ref::YAXArray, fut::YAXArray; detrend::Bool=false, P::Real=0.95, window::Int=15, rankn::Int=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, wet_threshold::Real=1.0, interp=Linear(), extrap=Interpolations.Flat(), gevparams::DataFrame=DataFrame(), frac=0.60, power=3.0, runlength::Int=2)
    obs_tdim = _time_dim_symbol(obs)
    ref_tdim = _time_dim_symbol(ref)
    fut_tdim = _time_dim_symbol(fut)

    qqmap_base = qqmap(obs, ref, fut; method="multiplicative", detrend=detrend, window=window, rankn=rankn, qmin=qmin, qmax=qmax, thresnan=thresnan, keep_original=keep_original, wet_threshold=wet_threshold, interp=interp, extrap=extrap)

    obs_noleap = drop29thfeb(obs, dimt=obs_tdim)
    ref_noleap = drop29thfeb(ref, dimt=ref_tdim)
    fut_noleap = drop29thfeb(fut, dimt=fut_tdim)

    datevec_ref = lookup(ref_noleap, ref_tdim)
    datevec_fut = lookup(fut_noleap, fut_tdim)
    movingwindow = (Base.maximum(year.(datevec_fut)) - Base.minimum(year.(datevec_fut))) > (Base.maximum(year.(datevec_ref)) - Base.minimum(year.(datevec_ref))) * 1.5

    spatial_dims = Symbol[d for d in name.(fut_noleap.axes) if d != fut_tdim]

    obs2d, _ = _cube_to_timeseries_matrix(obs_noleap, obs_tdim, spatial_dims)
    ref2d, _ = _cube_to_timeseries_matrix(ref_noleap, ref_tdim, spatial_dims)
    fut2d, _ = _cube_to_timeseries_matrix(fut_noleap, fut_tdim, spatial_dims)
    qq2d, qq_spatial_shape = _cube_to_timeseries_matrix(qqmap_base, :time, spatial_dims)
    dataout2d = copy(qq2d)

    latflat = Float64[]
    lonflat = Float64[]
    if !isempty(gevparams)
        latdim = _find_dim_in_axes(fut_noleap, (:lat, :latitude, :rlat, :y))
        londim = _find_dim_in_axes(fut_noleap, (:lon, :longitude, :rlon, :x))
        (latdim === nothing || londim === nothing) && error("`gevparams` was provided but latitude/longitude dimensions could not be inferred.")
        latflat = _flatten_coord_values(fut_noleap, spatial_dims, latdim)
        lonflat = _flatten_coord_values(fut_noleap, spatial_dims, londim)
    end

    for k in axes(obs2d, 1)
        obsvec = Float64.(coalesce.(obs2d[k, :], NaN))
        refvec = Float64.(coalesce.(ref2d[k, :], NaN))
        futvec = Float64.(coalesce.(fut2d[k, :], NaN))

        obs_nanfrac = count(isnan, obsvec) / length(obsvec)
        ref_nanfrac = count(isnan, refvec) / length(refvec)
        fut_nanfrac = count(isnan, futvec) / length(futvec)

        if obs_nanfrac < thresnan && ref_nanfrac < thresnan && fut_nanfrac < thresnan
            thres = _threshold_from_vectors(obsvec, refvec; P=P)
            !isfinite(thres) && continue

            gpd_obs = nothing
            if isempty(gevparams)
                gpd_obs, _, _ = _fit_gpd_distribution(obsvec, thres; runlength=runlength)
            else
                gpd_obs = estimate_gpd(latflat[k], lonflat[k], gevparams, thres)
            end
            gpd_obs === nothing && continue

            if movingwindow
                idxi, idxf, idxi_corr, idxf_corr = build_idx_biascorrect(refvec, futvec)

                for c in eachindex(idxi)
                    futvec_tmp = futvec[idxi[c]:idxf[c]]
                    gpd_fut, ex_idx, maxvalues = _fit_gpd_distribution(futvec_tmp, thres; runlength=runlength)
                    (gpd_fut === nothing || isempty(ex_idx)) && continue

                    fut_cdf = Extremes.cdf.(Ref(gpd_fut), maxvalues .- thres)
                    fut_cdf = clamp.(fut_cdf, 1e-6, 1 - 1e-6)
                    newfut = Extremes.quantile.(Ref(gpd_obs), fut_cdf) .+ thres

                    transition = _extremes_weight(futvec_tmp[ex_idx]; frac=frac, power=power)
                    ex_global = idxi[c] .+ ex_idx .- 1

                    in_corr = (ex_global .< idxf_corr[c]) .& (ex_global .> idxi_corr[c])
                    !any(in_corr) && continue

                    selected_global = ex_global[in_corr]
                    selected_new = newfut[in_corr]
                    selected_transition = transition[in_corr]
                    qqbase = qq2d[k, selected_global]

                    dataout2d[k, selected_global] .= selected_new .* selected_transition .+ qqbase .* (1 .- selected_transition)
                end
            else
                gpd_fut, ex_idx, maxvalues = _fit_gpd_distribution(futvec, thres; runlength=runlength)
                (gpd_fut === nothing || isempty(ex_idx)) && continue

                fut_cdf = Extremes.cdf.(Ref(gpd_fut), maxvalues .- thres)
                fut_cdf = clamp.(fut_cdf, 1e-6, 1 - 1e-6)
                newfut = Extremes.quantile.(Ref(gpd_obs), fut_cdf) .+ thres

                transition = _extremes_weight(futvec[ex_idx]; frac=frac, power=power)
                qqbase = qq2d[k, ex_idx]
                dataout2d[k, ex_idx] .= newfut .* transition .+ qqbase .* (1 .- transition)
            end
        end
    end

    target_names = collect(Symbol.(name.(qqmap_base.axes)))
    dataout = _matrix_to_cube_array(dataout2d, qq_spatial_shape, spatial_dims, :time, target_names)
    return YAXArray(qqmap_base.axes, dataout)
end


"""
    generate_random_perturbations(knots::AbstractVector, vec)

Generate random perturbations for a given vector based on the specified knots.

# Arguments
- `knots::AbstractVector`: The knot values used to group observations.
- `vec`: The values to perturb within each knot group.

# Returns
- A vector containing the perturbed values.

"""
function generate_random_perturbations(knots::AbstractVector, vec)
    @info "Generating perturbations..."
    # @info now()

    data = DataFrame(knots = knots, values = vec)

    observations = :values

    newdata = combine(DataFrames.groupby(data, :knots, sort = false), observations => add_perturbations => observations)

    return newdata[!,observations]
    
end

function add_perturbations(x)
    # Random vectors parameters
    init     = 0.000001
    interval = 0.00001
    stop     = 0.001
    R = init:interval:stop
    
    return x .+ rand(R, length(x))
end


"""
    get_seasons(months)

Given a list of months, returns a list of corresponding seasons.

# Arguments
- `months::Array{Int}`: A list of months represented as integers.

# Returns
- `seasons::Array{String}`: A list of corresponding seasons.


"""
function get_seasons(months)
    seasons = []
    for month in months
        if month in [12, 1, 2]
            push!(seasons, "Winter")
        elseif month in [3, 4, 5]
            push!(seasons, "Spring")
        elseif month in [6, 7, 8]
            push!(seasons, "Summer")
        elseif month in [9, 10, 11]
            push!(seasons, "Autumn")
        end
    end
    return seasons
end


"""
    build_itp(refP, sf_refP, interp, extrap)

Build an interpolation function using the reference quantiles and correction factors.

# Arguments
- `refP`: Array of reference quantiles.
- `sf_refP`: Correction factors or shifts associated with `refP`.
- `interp`: Interpolation method to use.
- `extrap`: Extrapolation method to use.

# Returns
- `itp`: Interpolation function.
"""
function build_itp(refP, sf_refP, interp, extrap)
    
    if length(unique(refP)) < length(refP)
        Interpolations.deduplicate_knots!(refP, move_knots = false)
    end
    itp = interpolate((refP,), sf_refP, Gridded(interp))          
    itp = extrapolate(itp, extrap) # add extrapolation
    
    return itp
    
end

"""
    drop29thfeb(ds; dimx=:lon, dimy=:lat, dimt=:time)

Removes February 29th from the given YAXArray. This function is used for bias correction.
    
# Arguments
- `ds`: The input YAXArray from which February 29th needs to be removed.
- `dimt`: The time dimension symbol.

# Returns
- `YAXArray`: The modified YAXArray with February 29th removed and recast on a no-leap time axis.
"""
function drop29thfeb(ds; dimx=:lon, dimy=:lat, dimt=:time)
    
    ds_subset = ds[Dim{dimt}(Where(x -> Dates.monthday(x) != (2,29)))]

    # Old time vector
    date_vec = lookup(ds_subset, dimt)
    # New time vector
    datevec_noleap = CFTime.reinterpret.(DateTimeNoLeap, date_vec)

    # Rebuild the YAXArray
    newarray = set(ds_subset, dimt => datevec_noleap)

    return newarray

    
end

"""
    find_julianday_idx(julnb, ijulian, window)

Find the indices of julnb that fall within a specified window around ijulian.

The selection wraps around the year boundary so that early-January and late-December days
can be grouped in the same moving window.

# Arguments
- `julnb::Array{Int,1}`: Array of julian day numbers.
- `ijulian::Int`: Julian day number around which to find indices.
- `window::Int`: Size of the window around ijulian.

# Returns
- `idx::BitArray{1}`: Boolean array indicating the indices that fall within the window.
"""
function find_julianday_idx(julnb, ijulian, window)
    if ijulian <= window
        idx = @. (julnb >= 1) & (julnb <= (ijulian + window)) | (julnb >= (365 - (window - ijulian))) & (julnb <= 365)
    elseif ijulian > 365 - window
        idx = @. (julnb >= 1) & (julnb <= ((window-(365-ijulian)))) | (julnb >= (ijulian-window)) & (julnb <= 365)
    else
        idx = @. (julnb <= (ijulian + window)) & (julnb >= (ijulian - window))
    end
    return idx
end

"""
    add_attribs(ds, ivar, METHODS, imethod)

Legacy helper that rebuilds a `YAXArray` after bias correction and selects basic
metadata handling based on the variable name.

# Arguments
- `ds`: The array to wrap as a `YAXArray`.
- `ivar`: Variable identifier, typically `"tasmin"`, `"tasmax"`, or `"pr"`.
- `METHODS`: Unused legacy argument kept for compatibility.
- `imethod`: Unused legacy argument kept for compatibility.

# Returns
- A `YAXArray` rebuilt from `ds`.

# Notes
- This helper currently preserves the input axes and does not attach variable-specific metadata beyond rebuilding the array.
"""
function add_attribs(ds, ivar, METHODS, imethod)
    
    if ivar == "tasmin"
        ds = YAXArray(ds.axes, ds)#, Dict(["standard_name" => "air_temperature", "units" => "degC"]))
        
    elseif ivar == "tasmax"
        ds = YAXArray(ds.axes, ds)#, Dict(["standard_name" => "air_temperature", "units" => "degC"]))
        
    elseif ivar == "pr"
        ds = YAXArray(ds.axes, ds)#, Dict(["standard_name" => "precipitation_amount", "units" => "mm/d"]))
    else
        @error "Wrong variable?!?"
    end

    return ds
end

"""
    polyfit(vec; order=4)

Fit a polynomial trend to `vec` using time indices `1:length(vec)`.

# Arguments
- `vec`: The input vector containing the data points to fit.
- `order`: The order of the polynomial to fit (default is 4).

# Returns
- A polynomial object fitted to `vec`.
"""
function polyfit(vec; order=4)

    x = 1:length(vec)
    
    polynomial = Polynomials.fit(x, vec, order)
    # polynomial[0] = 0.0 # pkoi 0.0?
    return polynomial
    
end

function _detrend_series(vec; order=4)
    x = collect(1:length(vec))
    finite_idx = findall(isfinite, vec)

    if isempty(finite_idx)
        return copy(vec), zeros(Float64, length(vec))
    end

    fit_order = min(order, length(finite_idx) - 1)
    polynomial = Polynomials.fit(x[finite_idx], vec[finite_idx], fit_order)
    trend = polynomial.(x)

    return vec .- trend, trend
end

# """
#     polyval(vec, polynomial)

#     Returns a vector containing the values, as estimated from polynomial function polynomial.
# """
# function polyval(vec, polynomial)

#     return polynomial.(vec) # fonction n'est plus nécessaire
    
# end
