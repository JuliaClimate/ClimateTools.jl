"""
    qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method="Additive", detrend=true, window::Int=15, rankn::Int=50, thresnan::Float64=0.1, keep_original::Bool=false, interp::Function = Linear(), extrap::Function = Flat())

Quantile-Quantile mapping bias correction. For each julian day of the year (+/- **window** size), a transfer function is estimated through an empirical quantile-quantile mapping.

The quantile-quantile transfer function between **ref** and **obs** is etimated on a julian day (and grid-point) basis with a moving window around the julian day. Hence, for a given julian day, the transfer function is then applied to the **fut** dataset for a given julian day.

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
function qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method::String="Additive", detrend::Bool=true, window::Int=15, rankn::Int=50, thresnan::Float64=0.1, keep_original::Bool=false, interp=Linear(), extrap=Flat())

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
    dataout = fill(convert(typeof(fut[1].data[1]), NaN), (size(fut[1], 1), size(fut[1],2), size(fut[1], 3)))::Array{typeof(fut[1].data[1]), T} where T

    if minimum(ref_jul) == 1 && maximum(ref_jul) == 365
        days = 1:365
    else
        days = minimum(ref_jul)+window:maximum(ref_jul)-window
        start = Dates.monthday(minimum(datevec_obs))
        finish = Dates.monthday(maximum(datevec_obs))
    end

    # Reshape. Allows multi-threading on grid points.
    obsin = reshape(obs[1].data, (size(obs[1].data, 1)*size(obs[1].data, 2), size(obs[1].data, 3)))
    refin = reshape(ref[1].data, (size(ref[1].data, 1)*size(ref[1].data, 2), size(ref[1].data, 3)))
    futin = reshape(fut[1].data, (size(fut[1].data, 1)*size(fut[1].data, 2), size(fut[1].data, 3)))
    dataoutin = reshape(dataout, (size(dataout, 1)*size(dataout, 2), size(dataout, 3)))

    # Looping over grid points using multiple-dispatch calls to qqmap
    Threads.@threads for k = 1:size(obsin, 1)

        obsvec = obsin[k,:]
        refvec = refin[k,:]
        futvec = futin[k,:]

        dataoutin[k, :] = qqmap(obsvec, refvec, futvec, days, obs_jul, ref_jul, fut_jul, method=method, detrend=detrend, window=window, rankn=rankn, thresnan=thresnan, keep_original=keep_original, interp=interp, extrap=extrap)

    end

    # Apply mask
    dataout = applymask(dataout, obs.msk)

    # Rebuilding ClimGrid
    lonsymbol = Symbol(fut.dimension_dict["lon"])
    latsymbol = Symbol(fut.dimension_dict["lat"])

    dataout2 = AxisArray(dataout, Axis{lonsymbol}(fut[1][Axis{lonsymbol}][:]), Axis{latsymbol}(fut[1][Axis{latsymbol}][:]),Axis{:time}(datevec_fut))

    C = ClimGrid(dataout2; longrid=fut.longrid, latgrid=fut.latgrid, msk=fut.msk, grid_mapping=fut.grid_mapping, dimension_dict=fut.dimension_dict, timeattrib=fut.timeattrib, model=fut.model, frequency=fut.frequency, experiment=fut.experiment, run=fut.run, project=fut.project, institute=fut.institute, filename=fut.filename, dataunits=fut.dataunits, latunits=fut.latunits, lonunits=fut.lonunits, variable=fut.variable, typeofvar=fut.typeofvar, typeofcal=fut.typeofcal, varattribs=fut.varattribs, globalattribs=fut.globalattribs)

    if detrend == true
        C = C + poly_values
    end

    return C

end

"""
    qqmap(obs::Array{N, 1} where N, ref::Array{N, 1} where N, fut::Array{N, 1} where N; method="Additive", detrend=true, window=15, rankn=50, thresnan=0.1, keep_original=false, interp::Function=Linear(), extrap::Function=Flat())

Quantile-Quantile mapping bias correction for single vector. This is a low level function used by qqmap(A::ClimGrid ..), but can work independently.

"""
function qqmap(obsvec::Array{N,1} where N, refvec::Array{N,1} where N, futvec::Array{N,1} where N, days, obs_jul, ref_jul, fut_jul; method::String="Additive", detrend::Bool=true, window::Int64=15, rankn::Int64=50, thresnan::Float64=0.1, keep_original::Bool=false, interp=Linear(), extrap=Flat())

    # range over which quantiles are estimated
    P = range(0.01, stop=0.99, length=rankn)
    obsP = similar(obsvec, length(P))
    refP = similar(refvec, length(P))
    sf_refP = similar(refvec, length(P))

    # Prepare output array
    dataout = similar(futvec, (size(futvec)))

    # LOOP OVER ALL DAYS OF THE YEAR
    # TODO Define "days" instead of 1:365
    for ijulian in days

        # idx for values we want to correct
        idxfut = (fut_jul .== ijulian)

        # Find all index of moving window around ijulian day of year
        idxobs = ClimateTools.find_julianday_idx(obs_jul, ijulian, window)
        idxref = ClimateTools.find_julianday_idx(ref_jul, ijulian, window)

        obsval = obsvec[idxobs]# .+ eps(1.0) # values to use as ground truth
        refval = refvec[idxref]# .+ eps(1.0)# values to use as reference for sim
        futval = futvec[idxfut]# .+ eps(1.0) # values to correct

        obsvalnan = isnan.(obsval)
        refvalnan = isnan.(refval)
        futvalnan = isnan.(futval)

        if (sum(obsvalnan) < (length(obsval) * thresnan)) & (sum(refvalnan) < (length(refval) * thresnan)) & (sum(futvalnan) < (length(futval) * thresnan))

            # Estimate quantiles for obs and ref for ijulian
            quantile!(obsP, obsval[.!obsvalnan], P, sorted=false)
            quantile!(refP, refval[.!refvalnan], P, sorted=false)

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
                dataout[idxfut] .= NaN
                # DO NOTHING (e.g. if there is no reference, we want NaNs and not original values)
            end
        end
    end

    return dataout

end
