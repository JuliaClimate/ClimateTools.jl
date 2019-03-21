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
    # Remove trend if specified

    # Consistency checks # TODO add more checks for grid definition
    @argcheck size(obs[1], 1) == size(ref[1], 1) == size(fut[1], 1)
    @argcheck size(obs[1], 2) == size(ref[1], 2) == size(fut[1], 2)


    if detrend == true
        # Obs
        obs = ClimateTools.dropfeb29(obs) # Removes 29th February
        obs_polynomials = ClimateTools.polyfit(obs)
        obs = obs - ClimateTools.polyval(obs, obs_polynomials)
        # Ref
        ref = ClimateTools.dropfeb29(ref) # Removes 29th February
        ref_polynomials = ClimateTools.polyfit(ref)
        ref = ref - ClimateTools.polyval(ref, ref_polynomials)
        # Fut
        fut = ClimateTools.dropfeb29(fut) # Removes 29th February
        fut_polynomials = ClimateTools.polyfit(fut)
        poly_values = ClimateTools.polyval(fut, fut_polynomials)
        fut = fut - poly_values
    end

    #Get date vectors
    datevec_obs = get_timevec(obs) # [1][Axis{:time}][:]
    datevec_ref = get_timevec(ref) # [1][Axis{:time}][:]
    datevec_fut = get_timevec(fut) # [1][Axis{:time}][:]

    # Modify dates (e.g. 29th feb are dropped/lost by default)
    # TODO simplify corrjuliandays. Unneccessary calls
    ~, ~, datevec_obs2 = ClimateTools.corrjuliandays(obs[1][1,1,:].data, datevec_obs)
    ~, ref_jul, ~ = ClimateTools.corrjuliandays(ref[1][1,1,:].data, datevec_ref)
    futvec2, ~, datevec_fut2 = ClimateTools.corrjuliandays(fut[1][1,1,:].data, datevec_fut)

    # Prepare output array (replicating data type of fut ClimGrid)
    dataout = fill(convert(typeof(fut[1].data[1]), NaN), (size(fut[1], 1), size(fut[1],2), size(futvec2, 1)))::Array{typeof(fut[1].data[1]), T} where T
    # dataout = fill(convert(typeof(fut[1].data[1]), NaN), (size(fut[1], 1), size(fut[1],2), size(futvec2, 1)))::Array{typeof(fut[1].data[1]), T} where T

    if minimum(ref_jul) == 1 && maximum(ref_jul) == 365
        days = 1:365
    else
        days = minimum(ref_jul)+window:maximum(ref_jul)-window
        start = Dates.monthday(minimum(datevec_obs2))
        finish = Dates.monthday(maximum(datevec_obs2))
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

        dataoutin[k, :] = qqmap(obsvec, refvec, futvec, days, datevec_obs, datevec_ref, datevec_fut, method=method, detrend=detrend, window=window, rankn=rankn, thresnan=thresnan, keep_original=keep_original, interp=interp, extrap=extrap)

    end

    lonsymbol = Symbol(fut.dimension_dict["lon"])
    latsymbol = Symbol(fut.dimension_dict["lat"])

    # Apply mask
    dataout = applymask(dataout, obs.msk)

    dataout2 = AxisArray(dataout, Axis{lonsymbol}(fut[1][Axis{lonsymbol}][:]), Axis{latsymbol}(fut[1][Axis{latsymbol}][:]),Axis{:time}(datevec_fut2))

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
function qqmap(obsvec::Array{N,1} where N, refvec::Array{N,1} where N, futvec::Array{N,1} where N, days, datevec_obs, datevec_ref, datevec_fut; method::String="Additive", detrend::Bool=true, window::Int64=15, rankn::Int64=50, thresnan::Float64=0.1, keep_original::Bool=false, interp=Linear(), extrap=Flat())

    # range over which quantiles are estimated
    P = range(0.01, stop=0.99, length=rankn)
    obsP = similar(obsvec, length(P))
    refP = similar(refvec, length(P))
    sf_refP = similar(refvec, length(P))


    # Get correct julian days (e.g. we can't have a mismatch of calendars between observed and models ref/fut. For instance, julian day 200 is not the same for a Standard calendar and a NoLeap calendar in a leap year)
    obsvec2, obs_jul, ~ = ClimateTools.corrjuliandays(obsvec, datevec_obs)
    refvec2, ref_jul, ~ = ClimateTools.corrjuliandays(refvec, datevec_ref)
    futvec2, fut_jul, ~ = ClimateTools.corrjuliandays(futvec, datevec_fut)

    # Prepare output array
    dataout = similar(futvec2, (size(futvec2)))

    # LOOP OVER ALL DAYS OF THE YEAR
    # TODO Define "days" instead of 1:365
    for ijulian in days

        # idx for values we want to correct
        idxfut = (fut_jul .== ijulian)

        # Find all index of moving window around ijulian day of year
        idxobs = ClimateTools.find_julianday_idx(obs_jul, ijulian, window)
        idxref = ClimateTools.find_julianday_idx(ref_jul, ijulian, window)

        # rng = MersenneTwister(1234)

        obsval = obsvec2[idxobs]# .+ eps(1.0) # values to use as ground truth
        refval = refvec2[idxref]# .+ eps(1.0)# values to use as reference for sim
        futval = futvec2[idxfut]# .+ eps(1.0) # values to correct

        if (sum(isnan.(obsval)) < (length(obsval) * thresnan)) & (sum(isnan.(refval)) < (length(refval) * thresnan)) & (sum(isnan.(futval)) < (length(futval) * thresnan))

            # Estimate quantiles for obs and ref for ijulian
            quantile!(obsP, obsval[.!isnan.(obsval)], P, sorted=false)
            quantile!(refP, refval[.!isnan.(refval)], P, sorted=false)

            if lowercase(method) == "additive" # used for temperature
                sf_refP .= obsP .- refP
                itp = interpolate((refP,), sf_refP, Gridded(interp))
                itp = extrapolate(itp, extrap) # add extrapolation
                # dataout[idxfut] .= itp(futval) .+ futval
                futval .= itp(futval) .+ futval

            elseif lowercase(method) == "multiplicative" # used for precipitation
                sf_refP .= obsP ./ refP
                sf_refP[sf_refP .< 0] .= eps(1.0)
                sf_refP[isnan.(sf_refP)] .= 1.0#eps(1.0)#1.0
                itp = interpolate((refP,), sf_refP, Gridded(interp))
                itp = extrapolate(itp, extrap) # add extrapolation
                futval .= itp(futval) .* futval
                # dataout[idxfut] .= itp(futval) .* futval

                futval[isnan.(futval)] .= 0.0

            else
                error("Wrong method")
            end
            # Replace values with new ones
            dataout[idxfut] .= futval
        else

            if keep_original
                # # Replace values with original ones (i.e. too may NaN values for robust quantile estimation)
                dataout[idxfut] = futval
            else
                dataout[idxfut] .= NaN
                # DO NOTHING (e.g. if there is no reference, we want NaNs and not original values)
            end
        end
    end

    # dataout[isnan.(dataout)] .= 0.0

    return dataout

end

# """
#     qqmaptf(obs::ClimGrid, ref::ClimGrid; partition::Float64 = 1.0, detrend::Bool=true, window::Int64=15, rankn::Int64=50, thresnan::Float64=0.1, keep_original::Bool=false, interp = Linear(), extrap = Flat())

# Transfer function based on quantile-quantile mapping bias correction. For each julian day, a transfer function is estimated through an empirical quantile-quantile mapping for the entire obs' ClimGrid extent. The quantile-quantile transfer function between **ref** and **obs** is etimated on a julian day basis with a moving window around the julian day. The transfer function can then be used to correct another dataset.

# **Options**
# partition::Float64 = 1.0. The proportion of grid-points (chosen randomly) used for the estimation of the transfer function. A transfer function is estimated for every chosen grid-points (and julian day) and averaged for the entire obs ClimGrid extent.

# **method::String = "Additive" (default) or "Multiplicative"**. Additive is used for most climate variables. Multiplicative is usually bounded variables such as precipitation and humidity.

# **detrend::Bool = true (default)**. A 4th order polynomial is adjusted to the time series and the residuals are corrected with the quantile-quantile mapping.

# **window::Int = 15 (default)**. The size of the window used to extract the statistical characteristics around a given julian day.

# **rankn::Int = 50 (default)**. The number of bins used for the quantile estimations. The quantiles uses by default 50 bins between 0.01 and 0.99. The bahavior between the bins is controlled by the interp keyword argument. The behaviour of the quantile-quantile estimation outside the 0.01 and 0.99 range is controlled by the extrap keyword argument.

# **interp = Interpolations.Linear() (default)**. When the data to be corrected lies between 2 quantile bins, the value of the transfer function is linearly interpolated between the 2 closest quantile estimation. The argument is from Interpolations.jl package.

# **extrap = Interpolations.Flat() (default)**. The bahavior of the quantile-quantile transfer function outside the 0.01-0.99 range. Setting it to Flat() ensures that there is no "inflation problem" with the bias correction. The argument is from Interpolation.jl package.

# """
# function qqmaptf(obs::ClimGrid, ref::ClimGrid; partition::Float64 = 1.0, method::String="Additive", detrend::Bool = true, window::Int64=15, rankn::Int64=50, interp = Linear(), extrap = Flat())
#     # Remove trend if specified
#     if detrend == true
#         obs = ClimateTools.dropfeb29(obs) # Removes 29th February
#         obs_polynomials = ClimateTools.polyfit(obs)
#         obs = obs - ClimateTools.polyval(obs, obs_polynomials)
#         ref = ClimateTools.dropfeb29(ref) # Removes 29th February
#         ref_polynomials = ClimateTools.polyfit(ref)
#         ref = ref - ClimateTools.polyval(ref, ref_polynomials)
#     end

#     # Checking if obs and ref are the same size
#     @argcheck size(obs[1], 1) == size(ref[1], 1)
#     @argcheck size(obs[1], 2) == size(ref[1], 2)

#     # range over which quantiles are estimated
#     P = range(0.01, stop=0.99, length=rankn)

#     #Get date vectors
#     datevec_obs = obs[1][Axis{:time}][:]
#     datevec_ref = ref[1][Axis{:time}][:]

#     # Modify dates (e.g. 29th feb are dropped/lost by default)
#     obsvec2, obs_jul, datevec_obs2 = ClimateTools.corrjuliandays(obs[1][1,1,:].data, datevec_obs)
#     refvec2, ref_jul, datevec_ref2 = ClimateTools.corrjuliandays(ref[1][1,1,:].data, datevec_ref)
#     if minimum(ref_jul) == 1 && maximum(ref_jul) == 365
#         days = 1:365
#     else
#         days = minimum(ref_jul)+window:maximum(ref_jul)-window
#         start = Dates.monthday(minimum(datevec_obs2))
#         finish = Dates.monthday(maximum(datevec_obs2))
#         warn(string("The reference ClimGrid doesn't cover all the year. The transfer function has been calculated from the ", minimum(ref_jul)+15, "th to the ", maximum(ref_jul)-15, "th julian day"))
#     end

#     # Number of points to sample
#     nx = round(Int, partition * size(obs[1], 1)) # Number of points in x coordinate
#     ny = round(Int, partition * size(obs[1], 2)) # Number of points in y coordinate
#     if nx == 0
#         nx = 1
#     end
#     if ny == 0
#         ny = 1
#     end
#     # Coordinates of the sampled points
#     x = sort(Random.randperm(size(obs[1],1))[1:nx])
#     y = sort(Random.randperm(size(obs[1],2))[1:ny])
#     # Make sure at least one point is not NaN
#     while isnan(obs[1][x[1],y[1],:].data[1])
#         x = sort(Random.randperm(size(obs[1],1))[1:nx])
#         y = sort(Random.randperm(size(obs[1],2))[1:ny])
#     end

#     # Create matrix of indices
#     X, Y = meshgrid(x, y)

#     # Initialization of the output
#     ITP = Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Interpolations.Gridded{typeof(interp)},Tuple{Array{Float64,1}},0},Interpolations.Gridded{typeof(interp)},Interpolations.OnGrid,typeof(extrap)}}(undef, 365)

#     # Loop over every julian days
#     println("Estimating transfer functions...This can take a while.")
#     # p = Progress(length(days), 1)
#     Threads.@threads for ijulian in days
#         # Index of ijulian ± window
#         idxobs = ClimateTools.find_julianday_idx(obs_jul, ijulian, window)
#         idxref = ClimateTools.find_julianday_idx(ref_jul, ijulian, window)
#         # Object containing observation/reference data of the n points on ijulian day
#         obsval = fill(NaN, sum(idxobs) * nx * ny)
#         refval = fill(NaN, sum(idxref) * nx * ny)

#         ipoint = 1
#         for (ix, iy) in zip(X, Y)
#             # for iy in y
#                 iobsvec2, iobs_jul, idatevec_obs2 = ClimateTools.corrjuliandays(obs[1][ix,iy,:].data, datevec_obs)
#                 irefvec2, iref_jul, idatevec_ref2 = ClimateTools.corrjuliandays(ref[1][ix,iy,:].data, datevec_ref)
#                 obsval[sum(idxobs)*(ipoint-1)+1:sum(idxobs)*ipoint] = iobsvec2[idxobs]
#                 refval[sum(idxref)*(ipoint-1)+1:sum(idxref)*ipoint] = irefvec2[idxref]
#                 ipoint += 1
#             # end
#         end

#         # Estimate quantiles for obs and ref for ijulian
#         obsP = quantile(obsval[.!isnan.(obsval)], P)
#         refP = quantile(refval[.!isnan.(refval)], P)
#         if lowercase(method) == "additive" # used for temperature
#             sf_refP = obsP - refP
#         elseif lowercase(method) == "multiplicative" # used for precipitation
#             sf_refP = obsP ./ refP
#             sf_refP[sf_refP .< 0] = 0.
#         end
#         # transfert function for ijulian
#         itp = interpolate((refP,), sf_refP, Gridded(interp))
#         itp = extrapolate(itp, extrap) # add extrapolation
#         ITP[ijulian] = itp
#         # next!(p)
#     end
#     ITPout = TransferFunction(ITP, method, detrend)
#     return ITPout
# end

# """
#     qqmap(fut::ClimGrid, ITP::TransferFunction)
#
# Quantile-Quantile mapping bias correction with a known transfer function. For each julian day of the year, use the right transfert function to correct *fut* values.
#
# """
# function qqmap(fut::ClimGrid, ITP::TransferFunction)
#     if ITP.detrend == true
#         fut = dropfeb29(fut) # Removes 29th February
#         fut_polynomials = polyfit(fut)
#         poly_values = polyval(fut, fut_polynomials)
#         fut = fut - poly_values
#     end
#     # Get date vectors
#     datevec_fut = fut[1][Axis{:time}][:]
#     futvec2, fut_jul, datevec_fut2 = corrjuliandays(fut[1][1,1,:].data, datevec_fut)
#     days = minimum(fut_jul):maximum(fut_jul)
#     # Prepare output array
#     dataout = fill(NaN, (size(fut[1], 1), size(fut[1],2), size(futvec2, 1)))::Array{N, T} where N where T
#
#     Threads.@threads for ijulian in days
#         idxfut = (fut_jul .== ijulian)
#         # Value to correct
#         # futval = futvec2[idxfut]
#         futval = fut[1][:,:,idxfut].data
#         # Transfert function for ijulian
#         itp = ITP.itp[ijulian]
#         # Correct futval
#         if lowercase(ITP.method) == "additive" # used for temperature
#             futnew = itp[futval] .+ futval
#         elseif lowercase(ITP.method) == "multiplicative" # used for precipitation
#             futnew = itp[futval] .* futval
#         else
#             error("Wrong method")
#         end
#         # futvec_corr[idxfut] = futnew
#         dataout[:,:,idxfut] = futnew
#     end
#
#     lonsymbol = Symbol(fut.dimension_dict["lon"])
#     latsymbol = Symbol(fut.dimension_dict["lat"])
#
#     dataout2 = AxisArray(dataout, Axis{lonsymbol}(fut[1][Axis{lonsymbol}][:]), Axis{latsymbol}(fut[1][Axis{latsymbol}][:]),Axis{:time}(datevec_fut2))
#
#     C = ClimGrid(dataout2; longrid=fut.longrid, latgrid=fut.latgrid, msk=fut.msk, grid_mapping=fut.grid_mapping, dimension_dict=fut.dimension_dict, timeattrib=fut.timeattrib, model=fut.model, frequency=fut.frequency, experiment=fut.experiment, run=fut.run, project=fut.project, institute=fut.institute, filename=fut.filename, dataunits=fut.dataunits, latunits=fut.latunits, lonunits=fut.lonunits, variable=fut.variable, typeofvar=fut.typeofvar, typeofcal=fut.typeofcal, varattribs=fut.varattribs, globalattribs=fut.globalattribs)
#
#     if ITP.detrend == true
#         C = C + poly_values
#     end
#
#     return C
# end
