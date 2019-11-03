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
function qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method::String = "Additive", detrend::Bool = true, window::Int = 15, rankn::Int = 50, thresnan::Float64 = 0.1, keep_original::Bool = false, interp = Linear(), extrap = Flat())

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

        obsvec = obsin[k,:]
        refvec = refin[k,:]
        futvec = futin[k,:]

        dataoutin[k, :] = qqmap(obsvec, refvec, futvec, days, obs_jul, ref_jul, fut_jul, method = method, detrend = detrend, window = window, rankn = rankn, thresnan = thresnan, keep_original = keep_original, interp = interp, extrap = extrap)

    end

    # Apply mask
    dataout = applymask(dataout, obs.msk)

    # Rebuilding ClimGrid
    lonsymbol = Symbol(fut.dimension_dict["lon"])
    latsymbol = Symbol(fut.dimension_dict["lat"])

    dataout2 = AxisArray(dataout, Axis{lonsymbol}(fut[1][Axis{lonsymbol}][:]), Axis{latsymbol}(fut[1][Axis{latsymbol}][:]), Axis{:time}(datevec_fut))

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
function qqmap(obsvec::Array{N,1} where N, refvec::Array{N,1} where N, futvec::Array{N,1} where N, days, obs_jul, ref_jul, fut_jul; method::String = "Additive", detrend::Bool = true, window::Int64 = 15, rankn::Int64 = 50, thresnan::Float64 = 0.1, keep_original::Bool = false, interp = Linear(), extrap = Flat())

    # range over which quantiles are estimated
    P = range(0.01, stop = 0.99, length = rankn)
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

        # obsval = obsvec[idxobs] .+ eps(1.0) # values to use as ground truth
        # refval = refvec[idxref] .+ eps(1.0)# values to use as reference for sim
        # futval = futvec[idxfut] .+ eps(1.0) # values to correct

        obsval = obsvec[idxobs] .+ rand(0.000001:0.00001:0.001, length(obsvec[idxobs])) # values to use as ground truth
        refval = refvec[idxref] .+ rand(0.000001:0.00001:0.001, length(refvec[idxref])) # values to use as reference for sim
        futval = futvec[idxfut] .+ rand(0.000001:0.00001:0.001, length(futvec[idxfut])) # values to correct

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
    biascorrect_extremes(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method::String="Multiplicative", detrend::Bool=true, window::Int=15, rankn::Int=50, thresnan::Float64=0.1, keep_original::Bool=false, interp=Linear(), extrap=Flat(), gev_params::DataFrame)

Correct the tail of the distribution with a paramatric method, using the parameters μ, σ and ξ contained in gevparams DataFrame.
"""
function biascorrect_extremes(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method::String = "Multiplicative", P::Real = 0.95, detrend::Bool = true, window::Int = 15, rankn::Int = 50, thresnan::Float64 = 0.1, keep_original::Bool = false, interp = Linear(), extrap = Flat(), gevparams::DataFrame)

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
    latgrid = reshape(obs.latgrid, (size(obs.latgrid, 1) * size(obs.latgrid, 2)))
    longrid = reshape(obs.longrid, (size(obs.longrid, 1) * size(obs.longrid, 2)))

    # Looping over grid points using multiple-dispatch calls to qqmap
    # Threads.@threads for k = 1:size(obsin, 1)
    for k = 1:size(obsin, 1)

        obsvec = obsin[k,:]
        refvec = refin[k,:]
        futvec = futin[k,:]

        # Estimate threshold
        thres = ClimateTools.get_threshold(obsvec, refvec, thres = P)

        # Find closest gev parameters
        l1 = (Float64(obs.latgrid[k]), Float64(obs.longrid[k]))
        l2 = (Array(gevparams[:lat]), Array(gevparams[:lon]))
        dist, idx = findmindist(l1, l2)

        # Extract GEV parameters
        μ = gevparams[idx[1], :mu]
        σ = gevparams[idx[1], :sigma]
        ξ = gevparams[idx[1], :xi]

        # Transform to GPD parameters and build reference GP distribution.
        μ₂ = ClimateTools.gev2gpd(μ, σ, ξ, thres)
        GPD_obs = Extremes.GeneralizedPareto(μ₂, σ, ξ)

        if maximum(year.(datevec_fut)) - minimum(year.(datevec_fut)) > (maximum(year.(datevec_ref)) - minimum(year.(datevec_ref)))*1.5
            # window_length will be roughly the length of REF
            ref_l = (maximum(year.(datevec_ref)) - minimum(year.(datevec_ref)) + 1)*365
            # in case of multiple members, we approximate ref length
            fut_l = length(futvec)
            w_length = fut_l/(floor(fut_l/ref_l))
            # The CON series is separated is a few time windows.
            # We also make the windows overlap to replicate a moving window
            R = w_length/2:w_length/3:fut_l
            for c in R
                idxi = round(Int, max(c-w_length/2+1,1))
                idxf = round(Int, min(c+w_length/2,fut_l))
                futvec_tmp = futvec[idxi:idxf]
                fut_jul_tmp = fut_jul[idxi:idxf]

                # However, we want each moving window to correct a smaller portion of the time series (PR. ??)
                idxi_corr = round(Int, w_length/2-(w_length/2)/3+1)
                if c == w_length/2
                    idxi_corr = 1
                end
                idxf_corr = round(Int, w_length/2+(w_length/2)/3)
                if c == maximum(w_length/2:w_length/3:fut_l) || idxi+idxf_corr-1 > fut_l
                    idxf_corr = min(round(Int, w_length-w_length/3),fut_l-idxi+1)
                end

                nbyears = round(Int,length(futvec_tmp)/365)
                N = nbyears*2 # By default, we correct roughly 2 values per year

                # refclusters = getcluster(refvec[k,:], thres, 1.0)
                futclusters = getcluster(futvec_tmp, thres, 1.0)

                # GPD_obs = gpdfit(obsclusters[:Max], threshold = thres)
                # GPD_ref = gpdfit(refclusters[:Max], threshold = thres)
                GPD_fut = gpdfit(futclusters[:Max], threshold = thres)

                # ref_cdf = Extremes.cdf.(GPD_ref, refclusters[:Max])
                fut_cdf = Extremes.cdf.(GPD_fut, futclusters[:Max])

                newfut = quantile.(GPD_obs, fut_cdf)

                datatmp = qqmap(obsvec, refvec, futvec_tmp, days, obs_jul, ref_jul, fut_jul_tmp, method=method, detrend=detrend, window=window, rankn=rankn, thresnan=thresnan, keep_original=keep_original, interp=interp, extrap=extrap)

                dataoutin[k, idxi_corr:idxf_corr] = datatmp[idxi_corr:idxf_corr]

                # TODO Linear interpolations between beginning of GPD and cdf = 0.75
                id_ex_corr = findall(futclusters[:Position] .< idxf_corr)

                dataoutin[k, futclusters[id_ex_corr,:Position]] = newfut[id_ex_corr]



            end
        else
            # refclusters = getcluster(refvec[k,:], thres, 1.0)
            futclusters = getcluster(futvec, thres, 1.0)
            GPD_fut = gpdfit(futclusters[:Max], threshold = thres)
            fut_cdf = Extremes.cdf.(GPD_fut, futclusters[:Max])
            newfut = quantile.(GPD_obs, fut_cdf)

            dataoutin[k, :] = qqmap(obsvec, refvec, futvec, days, obs_jul, ref_jul, fut_jul, method=method, detrend=detrend, window=window, rankn=rankn, thresnan=thresnan, keep_original=keep_original, interp=interp, extrap=extrap)


            dataoutin[k, futclusters[:Position]] = newfut


        end






        # TODO Add a moving window and estimate clusters and quantile accordingly
        # obsclusters = getcluster(obsin[k,:], thres)


    end

    # Apply mask
    dataout = applymask(dataout, obs.msk)

    # Rebuilding ClimGrid
    lonsymbol = Symbol(fut.dimension_dict["lon"])
    latsymbol = Symbol(fut.dimension_dict["lat"])

    dataout2 = AxisArray(dataout, Axis{lonsymbol}(fut[1][Axis{lonsymbol}][:]), Axis{latsymbol}(fut[1][Axis{latsymbol}][:]), Axis{:time}(datevec_fut))

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
