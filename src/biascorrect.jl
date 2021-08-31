"""
qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method="Additive", detrend=true, window::Int=15, rankn::Int=50, thresnan::Float64=0.1, keep_original::Bool=false, interp::Function = Linear(), extrap::Function = Flat())

Quantile-Quantile mapping bias correction.

For each julian day of the year (+/- **window** size), a transfer function is estimated through an empirical quantile-quantile mapping. The quantile-quantile transfer function between **ref** and **obs** is etimated on a julian day (and grid-point) basis with a moving window around the julian day. Hence, for a given julian day, the transfer function is then applied to the **fut** dataset for a given julian day.

**Options**

**method::String = "Additive" (default) or "Multiplicative"**. Additive is used for most climate variables. Multiplicative is usually bounded variables such as precipitation and humidity.

**detrend::Bool = true (default)**. A 4th order polynomial is adjusted to the time series and the residuals are corrected with the quantile-quantile mapping.

**window::Int = 15 (default)**. The size of the window used to extract the statistical characteristics around a given julian day.

**rankn::Int = 50 (default)**. The number of bins used for the quantile estimations. The quantiles uses by default 50 bins between 0.01 and 0.99. The bahavior between the bins is controlled by the interp keyword argument. The behaviour of the quantile-quantile estimation outside the 0.01 and 0.99 range is controlled by the extrap keyword argument.

**thresnan::Float64 = 0.1 (default)**. The fraction is missing values authorized for the estimation of the quantile-quantile mapping for a given julian days. If there is more than **treshnan** missing values, the output for this given julian days returns NaNs.

**keep_original::Bool = false (default)**. If **keep_original** is set to true, the values are set to the original values in presence of too many NaNs.

**interp = Interpolations.Linear() (default)**. When the data to be corrected lies between 2 quantile bins, the value of the transfer function is linearly interpolated between the 2 closest quantile estimation. The argument is from Interpolations.jl package.

**extrap = Interpolations.Flat() (default)**. The bahavior of the quantile-quantile transfer function outside the 0.01-0.99 range. Setting it to Flat() ensures that there is no "inflation problem" with the bias correction. The argument is from Interpolation.jl package.
"""
function qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method::String="Additive", detrend::Bool=true, window::Int=15, rankn::Int=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, interp=Linear(), extrap=Flat())

    # Consistency checks # TODO add more checks for grid definition
    @argcheck size(obs[1], 1) == size(ref[1], 1) == size(fut[1], 1)
    @argcheck size(obs[1], 2) == size(ref[1], 2) == size(fut[1], 2)

    # Modify dates (e.g. 29th feb are dropped/lost by default)
    obs = ClimateTools.dropfeb29(obs)
    ref = ClimateTools.dropfeb29(ref)
    fut = ClimateTools.dropfeb29(fut)

    # Remove trend if specified
    if detrend == true
        # Obs
        obs_polynomials = ClimateTools.polyfit(obs)
        obs = obs - ClimateTools.polyval(obs, obs_polynomials)
        # Ref
        ref_polynomials = ClimateTools.polyfit(ref)
        ref = ref - ClimateTools.polyval(ref, ref_polynomials)
        # Fut
        fut_polynomials = ClimateTools.polyfit(fut)
        poly_values = ClimateTools.polyval(fut, fut_polynomials)
        fut = fut - poly_values
    end

    #Get date vectors
    datevec_obs = get_timevec(obs)
    datevec_ref = get_timevec(ref)
    datevec_fut = get_timevec(fut)

    # Julian days vectors
    obs_jul = Dates.dayofyear.(datevec_obs)
    ref_jul = Dates.dayofyear.(datevec_ref)
    fut_jul = Dates.dayofyear.(datevec_fut)

    # Prepare output array (replicating data type of fut ClimGrid)
    dataout = fill(convert(typeof(fut[1].data[1]), NaN), (size(fut[1], 1), size(fut[1], 2), size(fut[1], 3)))::Array{typeof(fut[1].data[1]),T} where T

    if minimum(ref_jul) == 1 && maximum(ref_jul) == 365
        days = 1:365
    else
        days = minimum(ref_jul) + window:maximum(ref_jul) - window
        start = Dates.monthday(minimum(datevec_obs))
        finish = Dates.monthday(maximum(datevec_obs))
    end

    # Reshape. Allows multi-threading on grid points.
    obsin = reshape(obs[1].data, (size(obs[1].data, 1) * size(obs[1].data, 2), size(obs[1].data, 3)))
    refin = reshape(ref[1].data, (size(ref[1].data, 1) * size(ref[1].data, 2), size(ref[1].data, 3)))
    futin = reshape(fut[1].data, (size(fut[1].data, 1) * size(fut[1].data, 2), size(fut[1].data, 3)))
    dataoutin = reshape(dataout, (size(dataout, 1) * size(dataout, 2), size(dataout, 3)))

    # Looping over grid points using multiple-dispatch calls to qqmap
    Threads.@threads for k = 1:size(obsin, 1)
        # for k = 1:size(obsin, 1)

        obsvec = obsin[k,:]
        refvec = refin[k,:]
        futvec = futin[k,:]

        dataoutin[k, :] = qqmap(obsvec, refvec, futvec, days, obs_jul, ref_jul, fut_jul, method=method, detrend=detrend, window=window, rankn=rankn, qmin=qmin, qmax=qmax, thresnan=thresnan, keep_original=keep_original, interp=interp, extrap=extrap)

    end

    # Apply mask
    dataout = applymask(dataout, obs.msk)

    # Rebuilding ClimGrid
    lonsymbol = Symbol(fut.dimension_dict["lon"])
    latsymbol = Symbol(fut.dimension_dict["lat"])

    dataout2 = AxisArray(dataout, Axis{lonsymbol}(fut[1][Axis{lonsymbol}].val), Axis{latsymbol}(fut[1][Axis{latsymbol}].val), Axis{:time}(datevec_fut))

    timeattrib = fut.timeattrib
    timeattrib["calendar"] = "365_day"

    globalattribs = fut.globalattribs
    if haskey(globalattribs, "history")
        globalattribs["history"] = string(globalattribs["history"], " - ", "Bias-corrected with ClimateTools.jl")
    else
        globalattribs["history"] = string("Bias-corrected with ClimateTools.jl")
    end

    C = ClimGrid(dataout2; longrid = fut.longrid, latgrid = fut.latgrid, msk = fut.msk, grid_mapping = fut.grid_mapping, dimension_dict = fut.dimension_dict, timeattrib = timeattrib, model = fut.model, frequency = fut.frequency, experiment = fut.experiment, run = fut.run, project = fut.project, institute = fut.institute, filename = fut.filename, dataunits = fut.dataunits, latunits = fut.latunits, lonunits = fut.lonunits, variable = fut.variable, typeofvar = fut.typeofvar, typeofcal = "365_day", varattribs = fut.varattribs, globalattribs = globalattribs)

    if detrend == true
        C = C + poly_values
    end

    return C

end

"""
qqmap(obsvec::Array{N,1}, refvec::Array{N,1}, futvec::Array{N,1}, days, obs_jul, ref_jul, fut_jul; method::String="Additive", detrend::Bool=true, window::Int64=15, rankn::Int64=50, thresnan::Float64=0.1, keep_original::Bool=false, interp=Linear(), extrap=Flat())

Quantile-Quantile mapping bias correction for single vector. This is a low level function used by qqmap(A::ClimGrid ..), but can work independently.

"""
function qqmap(obsvec::Array{N,1} where N, refvec::Array{N,1} where N, futvec::Array{N,1} where N, days, obs_jul, ref_jul, fut_jul; method::String="Additive", detrend::Bool=true, window::Int64=15, rankn::Int64=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, interp=Linear(), extrap=Flat())

    # range over which quantiles are estimated
    P = range(qmin, stop=qmax, length = rankn)
    obsP = similar(obsvec, length(P))
    refP = similar(refvec, length(P))
    sf_refP = similar(refvec, length(P))

    # Prepare output array
    dataout = similar(futvec, (size(futvec)))

    # LOOP OVER ALL DAYS OF THE YEAR
    for ijulian in days

        # idx for values we want to correct
        idxfut = (fut_jul .== ijulian)

        # Find all index of moving window around ijulian day of year
        idxobs = ClimateTools.find_julianday_idx(obs_jul, ijulian, window)
        idxref = ClimateTools.find_julianday_idx(ref_jul, ijulian, window)

        obsval = obsvec[idxobs]
        refval = refvec[idxref]
        futval = futvec[idxfut]

        # Random perturbation is added to avoid division by zero
        if lowercase(method) == "multiplicative"
            # Seed for reproducible results
            Random.seed!(42)

            init     = 0.000001
            interval = 0.00001
            stop     = 0.001
            R = init:interval:stop

            obsval = obsval .+ rand(R, length(obsvec[idxobs]))
            refval = refval .+ rand(R, length(refvec[idxref]))
            futval = futval .+ rand(R, length(futvec[idxfut]))
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
                itp = interpolate((refP,), sf_refP, Gridded(interp))
                itp = extrapolate(itp, extrap) # add extrapolation
                # dataout[idxfut] .= itp(futval) .+ futval
                futval .= itp(futval) .+ futval

                futval[futvalnan] .= 0.0

            elseif lowercase(method) == "multiplicative" # used for precipitation
                sf_refP .= obsP ./ refP
                sf_refP[sf_refP .< 0] .= eps(1.0)
                sf_refP[isnan.(sf_refP)] .= 1.0#eps(1.0)#1.0
                itp = interpolate((refP,), sf_refP, Gridded(interp))
                itp = extrapolate(itp, extrap) # add extrapolation
                futval .= itp(futval) .* futval

                futval[futvalnan] .= 0.0

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

    return dataout

end

"""
    biascorrect_extremes(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; detrend=false, window::Int=15, rankn::Int=50, thresnan::Float64=0.1, keep_original::Bool=false, interp=Linear(), extrap=Flat(), gev_params::DataFrame, frac=0.25, power=1.0)

Combine the empirial Quantile-Quantile mapping (see [`qqmap`](@ref)) and Generalized Pareto Distribution bias-correction methods.

The tail of the distribution is corrected with a paramatric GPD, using provided parameters μ, σ and ξ contained in gevparams DataFrame. The DataFrame gevparams has column :lat, :lon, :mu, :sigma and :xi, representing the GEV parameters and the spatial location of the parameters. For each grid point, the function extracts the closest GEV parameters.

**Options specific to GPD**

**gevparams::DataFrame** represents the (external) GEV parameters to correct each grid-points.

**frac=0.25** is the fraction where the cutoff happens between QQM and GPD, as defined by **(maximum(x) - minimum(x))*frac** (e.g. For a maximum value of 150mm and a minimum value of 30mm, the linear transition will be between 30mm and 60mm).

**power=1.0** is the shape of the transition.

**Options specific to QQM**

**detrend::Bool = false (default)**. A 4th order polynomial is adjusted to the time series and the residuals are corrected.

**window::Int = 15 (default)**. The size of the window used to extract the statistical characteristics around a given julian day.

**rankn::Int = 50 (default)**. The number of bins used for the quantile estimations. The quantiles uses by default 50 bins between 0.01 and 0.99. The bahavior between the bins is controlled by the interp keyword argument. The behaviour of the quantile-quantile estimation outside the 0.01 and 0.99 range is controlled by the extrap keyword argument.

**thresnan::Float64 = 0.1 (default)**. The fraction is missing values authorized for the estimation of the quantile-quantile mapping for a given julian days. If there is more than **treshnan** missing values, the output for this given julian days returns NaNs.

**keep_original::Bool = false (default)**. If **keep_original** is set to true, the values are set to the original values in presence of too many NaNs.

**interp = Interpolations.Linear() (default)**. When the data to be corrected lies between 2 quantile bins, the value of the transfer function is linearly interpolated between the 2 closest quantile estimation. The argument is from Interpolations.jl package.

**extrap = Interpolations.Flat() (default)**. The bahavior of the quantile-quantile transfer function outside the 0.01-0.99 range. Setting it to Flat() ensures that there is no "inflation problem" with the bias correction. The argument is from Interpolation.jl package.

See also: [`qqmap`](@ref)
"""
function biascorrect_extremes(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; detrend=false, P::Real=0.95, window::Int=15, rankn::Int=50, qmin::Real=0.01, qmax::Real=0.99, thresnan::Float64=0.1, keep_original::Bool=false, interp=Linear(), extrap=Flat(), gevparams::DataFrame=DataFrame(), frac=0.25, power=1.0)

    # Consistency checks # TODO add more checks for grid definition
    @argcheck size(obs[1], 1) == size(ref[1], 1) == size(fut[1], 1)
    @argcheck size(obs[1], 2) == size(ref[1], 2) == size(fut[1], 2)

    qqmap_base = qqmap(obs, ref, fut, method="multiplicative", detrend=detrend, window=window, rankn=rankn, qmin=qmin, qmax=qmax, thresnan=thresnan, keep_original=keep_original, interp=interp, extrap=extrap)

    # Modify dates (e.g. 29th feb are dropped/lost by default)
    obs = ClimateTools.dropfeb29(obs)
    ref = ClimateTools.dropfeb29(ref)
    fut = ClimateTools.dropfeb29(fut)

    # Remove trend if specified
    if detrend == true
        # Obs
        obs_polynomials = ClimateTools.polyfit(obs)
        obs = obs - ClimateTools.polyval(obs, obs_polynomials)
        # Ref
        ref_polynomials = ClimateTools.polyfit(ref)
        ref = ref - ClimateTools.polyval(ref, ref_polynomials)
        # Fut
        fut_polynomials = ClimateTools.polyfit(fut)
        poly_values = ClimateTools.polyval(fut, fut_polynomials)
        fut = fut - poly_values
    end

    #Get date vectors
    datevec_obs = get_timevec(obs)
    datevec_ref = get_timevec(ref)
    datevec_fut = get_timevec(fut)

    # Condition for a moving window
    movingwindow = maximum(year.(datevec_fut)) - minimum(year.(datevec_fut)) > (maximum(year.(datevec_ref)) - minimum(year.(datevec_ref)))*1.5

    # Prepare output array (replicating data type of fut ClimGrid)
    # dataout = fill(convert(typeof(fut[1].data[1]), NaN), (size(fut[1], 1), size(fut[1], 2), size(fut[1], 3)))::Array{typeof(fut[1].data[1]),T} where T
    dataout = qqmap_base[1].data

    # Reshape. Allows multi-threading on grid points.
    obsin = reshape(obs[1].data, (size(obs[1].data, 1) * size(obs[1].data, 2), size(obs[1].data, 3)))
    refin = reshape(ref[1].data, (size(ref[1].data, 1) * size(ref[1].data, 2), size(ref[1].data, 3)))
    futin = reshape(fut[1].data, (size(fut[1].data, 1) * size(fut[1].data, 2), size(fut[1].data, 3)))
    dataoutin = reshape(dataout, (size(dataout, 1) * size(dataout, 2), size(dataout, 3)))
    latgrid = reshape(obs.latgrid, (size(obs.latgrid, 1) * size(obs.latgrid, 2)))
    longrid = reshape(obs.longrid, (size(obs.longrid, 1) * size(obs.longrid, 2)))
    qqvecin = reshape(qqmap_base[1].data, (size(qqmap_base[1].data, 1) * size(qqmap_base[1].data, 2), size(qqmap_base[1].data, 3)))

    # Looping over grid points using multiple-dispatch calls to qqmap
    # Threads.@threads for k = 1:size(obsin, 1)
    for k = 1:size(obsin, 1)

        obsvec = obsin[k,:]
        refvec = refin[k,:]
        futvec = futin[k,:]

        obsvalnan = isnan.(obsvec)
        refvalnan = isnan.(refvec)
        futvalnan = isnan.(futvec)

        if (sum(obsvalnan) < (length(obsvec) * thresnan)) & (sum(refvalnan) < (length(refvec) * thresnan)) & (sum(futvalnan) < (length(futvec) * thresnan))

            # Estimate threshold
            thres = ClimateTools.get_threshold(obsvec, refvec, thres=P)

            if isempty(gevparams)
                # If empty, we estimate the GPD on the reference dataset
                obsclusters = getcluster(obsvec, thres, 1.0)
                # Estimate GPD parameters of fut clusters
                GPD_obs = gpdfit(obsclusters[!,:Max], threshold = thres)

            else
                GPD_obs = ClimateTools.estimate_gpd(latgrid[k], longrid[k], gevparams, thres)
            end

            if movingwindow

                idxi, idxf, idxi_corr, idxf_corr = build_idx_biascorrect(refvec, futvec)

                for c = 1:length(idxi)

                    futvec_tmp = futvec[idxi[c]:idxf[c]]

                    # Get clusters
                    futclusters = getcluster(futvec_tmp, thres, 1.0)
                    # Estimate GPD parameters of fut clusters
                    GPD_fut = gpdfit(futclusters[!,:Max], threshold = thres)
                    # Calculate quantile values of maximum clusters values
                    fut_cdf = Extremes.cdf.(GPD_fut, futclusters[!,:Max])
                    # Estimation of fut values based on provided GPD parameters
                    newfut = quantile.(GPD_obs, fut_cdf)

                    # Linear transition over a quarter of the distance
                    exIDX = futclusters[!,:Position]

                    # Get the weights based on frac and power values
                    transition = ClimateTools.weight(futvec[exIDX], frac=frac, power=power)

                    # Apply linear transition
                    newfut_trans = (newfut .* transition) .+ (qqvecin[k, idxi[c] .+ exIDX .- 1] .* (1.0 .- transition))

                    # We put the new data in the bias-corrected vector
                    exIDX_corr = findall((idxi[c] .+ futclusters[!,:Position] .-1 .< idxf_corr[c]) .& (idxi[c] .+ futclusters[!,:Position] .- 1 .> idxi_corr[c]))
                    # exIDX_corr = findall(futclusters[!,:Position] .< idxf_corr)
                    # dataoutin[k, idxi+idxi_corr-1:idxi+idxf_corr-1] = qqvecin[k, idxi_corr:idxf_corr]
                    dataoutin[k, idxi[c] .+ futclusters[exIDX_corr,:Position] .- 1] .= newfut_trans[exIDX_corr]

                end
            else
                # refclusters = getcluster(refvec[k,:], thres, 1.0)
                futclusters = getcluster(futvec, thres, 1.0)
                GPD_fut = gpdfit(futclusters[!,:Max], threshold = thres)
                fut_cdf = Extremes.cdf.(GPD_fut, futclusters[!,:Max])
                newfut = quantile.(GPD_obs, fut_cdf)

                # Linear transition over a quarter of the distance
                exIDX = futclusters[!,:Position]

                # Get the weights based on frac and power values
                transition = ClimateTools.weight(futvec[exIDX], frac=frac, power=power)

                # Apply linear transition
                newfut_trans = (newfut .* transition) .+ (qqvecin[k, exIDX] .* (1.0 .- transition))

                # We put the new data in the bias-corrected vector
                # dataoutin[k, :] = qqvecin[k, :]
                dataoutin[k, futclusters[!,:Position]] = newfut_trans

            end

        # else
        #
        #     if keep_original
        #         # # Replace values with original ones (i.e. too may NaN values for robust quantile estimation)
        #         dataoutin[k, idxfut] .= futvec
        #     else
        #         # DO NOTHING (e.g. if there is no reference, we want NaNs and not original values)
        #         dataoutin[k, idxfut] .= NaN
        #     end
        end
    end

    # Apply mask
    dataout = applymask(dataout, obs.msk)

    # Rebuilding ClimGrid
    lonsymbol = Symbol(fut.dimension_dict["lon"])
    latsymbol = Symbol(fut.dimension_dict["lat"])

    dataout2 = AxisArray(dataout, Axis{lonsymbol}(fut[1][Axis{lonsymbol}].val), Axis{latsymbol}(fut[1][Axis{latsymbol}].val), Axis{:time}(datevec_fut))

    timeattrib = fut.timeattrib
    timeattrib["calendar"] = "365_day"

    globalattribs = fut.globalattribs
    if haskey(globalattribs, "history")
        globalattribs["history"] = string(globalattribs["history"], " - ", "Bias-corrected with ClimateTools.jl")
    else
        globalattribs["history"] = string("Bias-corrected for extreme values with ClimateTools.jl")
    end

    C = ClimGrid(dataout2; longrid = fut.longrid, latgrid = fut.latgrid, msk = fut.msk, grid_mapping = fut.grid_mapping, dimension_dict = fut.dimension_dict, timeattrib = timeattrib, model = fut.model, frequency = fut.frequency, experiment = fut.experiment, run = fut.run, project = fut.project, institute = fut.institute, filename = fut.filename, dataunits = fut.dataunits, latunits = fut.latunits, lonunits = fut.lonunits, variable = fut.variable, typeofvar = fut.typeofvar, typeofcal = "365_day", varattribs = fut.varattribs, globalattribs = globalattribs)

    if detrend == true
        C = C + poly_values
    end

    return C

end

"""
    gev2gpd(μ, σ, ξ, thres)

Transform GEV parameters to GPD parameters (Coles, 2001)
"""
function gev2gpd(μ, σ, ξ, thres)

    return σ₂ = σ + ξ * (thres - μ)
end


"""
    weight(x; frac=0.25, power=1.0)

Returns the associated weight used in weighting Quantile-Quantile mapping (QQM) and Generalized Pareto Distribution (GPD) bias-correction methods. Used internally by biascorrect_extremes.

*frac* is the fraction where the cutoff happens between QQM and GPD, as defined by **(maximum(x) - minimum(x))*frac** (e.g. For a maximum value of 150mm and a minimum value of 30mm, the linear transition will be between 30mm and 60mm) and power is the shape of the transition.
"""
function weight(x; frac=0.25, power=1.0)
    ω = ((x .- minimum(x))/((maximum(x) - minimum(x))*frac)).^power
    ω[ω .> 1.0] .= 1.0
    return ω
end


function estimate_gpd(lat, lon, gevparams, thres)

    # Find closest gev parameters
    l1 = (Float64(lat), Float64(lon))
    l2 = (Array(gevparams[!,:lat]), Array(gevparams[!,:lon]))
    dist, idx = findmindist(l1, l2)

    # Extract GEV parameters
    μ = gevparams[idx[1], :mu]
    σ = gevparams[idx[1], :sigma]
    ξ = gevparams[idx[1], :xi]

    # Transform to GPD parameters and build reference GP distribution.
    μ₂ = ClimateTools.gev2gpd(μ, σ, ξ, thres)
    GPD_obs = Extremes.GeneralizedPareto(μ₂, σ, ξ)

    return GPD_obs

end

function build_idx_biascorrect(refvec, futvec)

    # window_length will be roughly the length of REF
    w_length = length(refvec)#(maximum(year.(datevec_ref)) - minimum(year.(datevec_ref)) + 1)*365
    # in case of multiple members, we approximate ref length
    fut_l = length(futvec)

    # The fut series is separated is a few time windows.
    # We also make the windows overlap to replicate a moving window
    R = w_length/2:w_length/3:fut_l
    corrR = range(1, stop=fut_l, length=length(R)-1)

    idxi = Array{Int64}(undef, length(R))
    idxf = Array{Int64}(undef, length(R))
    idxi_corr = Array{Int64}(undef, length(R))
    idxf_corr = Array{Int64}(undef, length(R))

    for c = 1:length(R)

        # idxi = nodes_est[c,1]
        # idxf = nodes_est[c,2]
        # idxi_corr = nodes_corr[c]
        # idxf_corr = nodes_corr[c+1]

        idxi[c] = round(Int, max(R[c]-w_length/2+1,1))
        idxf[c] = round(Int, min(R[c]+w_length/2,fut_l))
        # futvec_tmp = futvec[idxi:idxf];

        # However, we want each moving window to correct a smaller portion of the time series
        idxi_corr[c] = round(Int, idxi[c]+(w_length/2-(w_length/2)/3) - 1)
        # idxi_corr[c] = round(Int, corrR[c])
        if c == 1
            idxi_corr[c] = 1
        end
        # idxf_corr = round(Int, idxi + (w_length/2+(w_length/2)/3) - 1)

        if c == length(R)#maximum(w_length/2:w_length/3:fut_l) || idxi+idxf_corr-1 > fut_l
            idxf_corr[c] = max(round(Int, idxi[c] + w_length-w_length/3),fut_l)
        else
            idxf_corr[c] = round(Int, idxi[c] + (w_length/2+(w_length/2)/3) - 1)
            # idxf_corr[c] = round(Int,corrR[c+1])
        end
        # @show idxi
        # @show idxf
        # @show idxi_corr
        # @show idxf_corr
    end

    idxi_corr[2:end] .-= 1

    return idxi, idxf, idxi_corr, idxf_corr
end




#     w_length = length(refvec)
#     fut_l = length(futvec)
#     nodes_corr = round.(Int,w_length/2:w_length/3:fut_l-w_length/2)
#     N = Array{Int64}(undef, length(nodes_corr)+2)
#     N[1] = 1
#     N[end] = fut_l
#     N[2:end-1] .= nodes_corr

#     N_est = Array{Int64}(undef, length(N), 2)

#     for inode = 1:length(N)

#         if inode == 1
#             N_est[inode, 1] = N[1]
#             N_est[inode, 2] = N[3]
#         elseif inode == length(N)
#             N_est[inode, 1] = N[end-1]
#             N_est[inode, 2] = N[end]
#         # elseif inode == length(N) - 1
#         #     N_est[inode, 1] = N[end-3]
#         #     N_est[inode, 2] = N[end-1]
#         else
#             N_est[inode, 1] = N[inode - 1]
#             N_est[inode, 2] = N[inode + 1]
#         end
#     end

#     return N, N_est
# end



#     # N[1] = 1
#     # N[end] = fut_l

#     for inode = 1:length(nodes_corr)
#         if inode == 1
#             N[inode, 1] = 1
#             N[inode, 2] = nodes_corr[inode]
#         elseif inode == length(nodes_corr)
#             N[inode] = fut_l
#         else
#             N[inode] = 0
#         end

#     end

#     nodes_corr = collect(Base.Iterators.partition(w_length/2:fut_l, round(Int,w_length/3)))
#     nodes_est = collect(Base.Iterators.partition(1:fut_l, round(Int,w_length)))

#     Ni = Vector{Int64}(undef, length(nodes_corr))
#     Nf = Vector{Int64}(undef, length(nodes_corr))

#     for inode = 1:length(nodes_corr)

#         for icorr = 1:length(nodes_est)

#             idxi = findall(nodes_corr[inode][1] .== nodes_est[icorr])
#             idxf = findall(nodes_corr[inode][end] .== nodes_est[icorr])

#             if !isempty(idxi)
#                 Ni[inode] = icorr
#             end
#             if !isempty(idxf)
#                 Nf[inode] = icorr
#             end
#         end
#     end

#     # @assert sum(Ni - Nf) == 0

#     return nodes_corr, nodes_est, Ni, Nf
# end


# for n = 1:30000
#     refvec = rand(10950)
#     futvec = rand(10950 + n)

#     nodes_corr, nodes_est, Ni, Nf = build_idx_biascorrect(refvec, futvec)
# end


# """
#     biascorrect_extremes(obs::AxisArray, ref::AxisArray, fut::AxisArray, μ, σ, ξ)
#
# Correct the tail of the distribution with a paramatric method, using the parameters μ, σ and ξ.
# """
# function biascorrect_extremes(obsvec::AxisArray, refvec::AxisArray, futvec::AxisArray, days, μ::Real, σ::Real, ξ::Real)
#
#     obs_jul = dayofyear.(obsvec[Axis{:time}][:])
#     ref_jul = dayofyear.(refvec[Axis{:time}][:])
#     fut_jul = dayofyear.(futvec[Axis{:time}][:])
#
#     # Prepare output array (replicating data type of fut ClimGrid)
#     dataout = fill(convert(typeof(futvec.data[1]), NaN), length(futvec))::Array{typeof(fut[1].data[1]),T} where T
#
#     if minimum(ref_jul) == 1 && maximum(ref_jul) == 365
#         days = 1:365
#     else
#         days = minimum(ref_jul) + window:maximum(ref_jul) - window
#         start = Dates.monthday(minimum(datevec_obs))
#         finish = Dates.monthday(maximum(datevec_obs))
#     end
#
#     threshold = mean([quantile(obsvec[obsvec .>= 1.0], P) quantile(refvec[refvec .>= 1.0], P)])
#
#     idx_obs = obsvec .> threshold
#     idx_ref = refvec .> threshold
#     idx_fut = futvec .> threshold
#
#     dataout .= qqmap(obsvec.data, refvec.data, futvec.data, days, obs_jul, ref_jul, fut_jul, method = "multiplicative", detrend = false)
#
#
# end
