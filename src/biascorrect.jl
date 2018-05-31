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

function qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method::String="Additive", detrend::Bool=true, window::Int=15, rankn::Int=50, thresnan::Float64=0.1, keep_original::Bool=false, interp = Linear(), extrap = Flat())

    # Consistency checks # TODO add more checks for grid definition
    @argcheck size(obs[1], 1) == size(ref[1], 1) == size(fut[1], 1)
    @argcheck size(obs[1], 2) == size(ref[1], 2) == size(fut[1], 2)

    #Get date vectors
    datevec_obs = obs[1][Axis{:time}][:]
    datevec_ref = ref[1][Axis{:time}][:]
    datevec_fut = fut[1][Axis{:time}][:]

    # Modify dates (e.g. 29th feb are dropped/lost by default)
    obsvec2, obs_jul, datevec_obs2 = corrjuliandays(obs[1][1,1,:].data, datevec_obs)
    refvec2, ref_jul, datevec_ref2 = corrjuliandays(ref[1][1,1,:].data, datevec_ref)
    futvec2, fut_jul, datevec_fut2 = corrjuliandays(fut[1][1,1,:].data, datevec_fut)
    # obsvec2, refvec2, futvec2, obs_jul, ref_jul, fut_jul, datevec_obs2, datevec_ref2, datevec_fut2 = corrjuliandays(obs[1][1,1,:].data, ref[1][1,1,:].data, fut[1][1,1,:].data, datevec_obs, datevec_ref, datevec_fut)

    # Prepare output array
    dataout = fill(NaN, (size(fut[1], 1), size(fut[1],2), size(futvec2, 1)))::Array{N, T} where N where T
    # dataout = fill(NaN, size(futvec2))::Array{N, T} where N where T

    p = Progress(size(obs[1], 3), 5)

    for k = 1:size(obs[1], 2)
        for j = 1:size(obs[1], 1)

            obsvec = obs[1][j,k,:].data
            refvec = ref[1][j,k,:].data
            futvec = fut[1][j,k,:].data

            dataout[j,k,:] = qqmap(obsvec, refvec, futvec, datevec_obs, datevec_ref, datevec_fut, method=method, detrend=detrend, window=window, rankn=rankn, thresnan=thresnan, keep_original=keep_original, interp=interp, extrap = extrap)

        end
        next!(p)
    end

    lonsymbol = Symbol(fut.dimension_dict["lon"])
    latsymbol = Symbol(fut.dimension_dict["lat"])

    dataout2 = AxisArray(dataout, Axis{lonsymbol}(fut[1][Axis{lonsymbol}][:]), Axis{latsymbol}(fut[1][Axis{latsymbol}][:]),Axis{:time}(datevec_fut2))

    return ClimGrid(dataout2, longrid=fut.longrid, latgrid=fut.latgrid, model=fut.model, project=fut.project, institute=fut.institute, experiment=fut.experiment, run=fut.run, filename=fut.filename, dataunits=fut.dataunits, latunits=fut.latunits, lonunits=fut.lonunits, variable=fut.variable, typeofvar=fut.typeofvar, typeofcal=fut.typeofcal, varattribs=fut.varattribs, globalattribs=fut.varattribs)

end


"""
    qqmap(obs::Array{N, 1} where N, ref::Array{N, 1} where N, fut::Array{N, 1} where N; method="Additive", detrend=true, window=15, rankn=50, thresnan=0.1, keep_original=false, interp::Function = Linear(), extrap::Function = Flat())

Quantile-Quantile mapping bias correction for single vector. This is a low level function used by qqmap(A::ClimGrid ..), but can work independently.

"""

function qqmap(obsvec::Array{N, 1} where N, refvec::Array{N, 1} where N, futvec::Array{N, 1} where N, datevec_obs, datevec_ref, datevec_fut; method::String="Additive", detrend::Bool=true, window::Int64=15, rankn::Int64=50, thresnan::Float64=0.1, keep_original::Bool=false, interp = Linear(), extrap = Flat())

    # range over which quantiles are estimated
    P = linspace(0.01, 0.99, rankn)

    # Get correct julian days (e.g. we can't have a mismatch of calendars between observed and models ref/fut)
    obsvec2, obs_jul, datevec_obs2 = corrjuliandays(obsvec, datevec_obs)
    refvec2, ref_jul, datevec_ref2 = corrjuliandays(refvec, datevec_ref)
    futvec2, fut_jul, datevec_fut2 = corrjuliandays(futvec, datevec_fut)
    # obsvec2, refvec2, futvec2, obs_jul, ref_jul, fut_jul, datevec_obs2, datevec_ref2, datevec_fut2 = corrjuliandays(obsvec, refvec, futvec, datevec_obs, datevec_ref, datevec_fut)

    # Prepare output array
    dataout = similar(futvec2, (size(futvec2)))

    # LOOP OVER ALL DAYS OF THE YEAR
    Threads.@threads for ijulian = 1:365

        # idx for values we want to correct
        idxfut = (fut_jul .== ijulian)

        # Find all index of moving window around ijulian day of year
        idxobs = find_julianday_idx(obs_jul, ijulian, window)
        idxref = find_julianday_idx(ref_jul, ijulian, window)

        obsval = obsvec2[idxobs] # values to use as ground truth
        refval = refvec2[idxref] # values to use as reference for sim
        futval = futvec2[idxfut] # values to correct

        if (sum(isnan.(obsval)) < (length(obsval) * thresnan)) & (sum(isnan.(refval)) < (length(refval) * thresnan)) & (sum(isnan.(futval)) < (length(futval) * thresnan))

            # Estimate quantiles for obs and ref for ijulian
            obsP = quantile(obsval[.!isnan.(obsval)], P)
            refP = quantile(refval[.!isnan.(refval)], P)

            if lowercase(method) == "additive" # used for temperature
                sf_refP = obsP - refP
                itp = interpolate((refP,), sf_refP, Gridded(interp))
                itp = extrapolate(itp, extrap) # add extrapolation
                futnew = itp[futval] + futval
            elseif lowercase(method) == "multiplicative" # used for precipitation
                sf_refP = obsP ./ refP
                sf_refP[sf_refP .< 0] = 0.
                itp = interpolate((refP,), sf_refP, Gridded(interp))
                itp = extrapolate(itp, extrap) # add extrapolation
                futnew = itp[futval] .* futval
            else
                error("Wrong method")
            end
            # Replace values with new ones
            dataout[idxfut] = futnew
        else

            if keep_original
                # # Replace values with original ones (i.e. too may NaN values for robust quantile estimation)
                dataout[idxfut] = futval
            else
                dataout[idxfut] = NaN
                # DO NOTHING (e.g. if there is no reference, we want NaNs and not original values)
            end
        end
    end

    return dataout

end

"""
    qqmaptf(obs::ClimGrid, ref::ClimGrid; partition::Float64 = 1.0, detrend::Bool=true, window::Int64=15, rankn::Int64=50, thresnan::Float64=0.1, keep_original::Bool=false, interp = Linear(), extrap = Flat())

Transfer function based on quantile-quantile mapping bias correction. For each julian day, a transfer function is estimated through an empirical quantile-quantile mapping for the entire obs' ClimGrid extent. The quantile-quantile transfer function between **ref** and **obs** is etimated on a julian day basis with a moving window around the julian day. The transfer function can then be used to correct another dataset.

**Options**
partition::Float64 = 1.0. The proportion of grid-points (chosen randomly) used for the estimation of the transfer function. A transfer function is estimated for every chosen grid-points (and julian day) and averaged for the entire obs ClimGrid extent.

**method::String = "Additive" (default) or "Multiplicative"**. Additive is used for most climate variables. Multiplicative is usually bounded variables such as precipitation and humidity.

**window::Int = 15 (default)**. The size of the window used to extract the statistical characteristics around a given julian day.

**rankn::Int = 50 (default)**. The number of bins used for the quantile estimations. The quantiles uses by default 50 bins between 0.01 and 0.99. The bahavior between the bins is controlled by the interp keyword argument. The behaviour of the quantile-quantile estimation outside the 0.01 and 0.99 range is controlled by the extrap keyword argument.

**interp = Interpolations.Linear() (default)**. When the data to be corrected lies between 2 quantile bins, the value of the transfer function is linearly interpolated between the 2 closest quantile estimation. The argument is from Interpolations.jl package.

**extrap = Interpolations.Flat() (default)**. The bahavior of the quantile-quantile transfer function outside the 0.01-0.99 range. Setting it to Flat() ensures that there is no "inflation problem" with the bias correction. The argument is from Interpolation.jl package.
"""
# TODO what happen when there is a lot of NaNs.
function qqmaptf(obs::ClimGrid, ref::ClimGrid; partition::Float64 = 1.0, method::String="Additive", window::Int64=15, rankn::Int64=50, interp = Linear(), extrap = Flat())
    # Checking if obs and ref are the same size
    @argcheck size(obs[1], 1) == size(ref[1], 1)
    @argcheck size(obs[1], 2) == size(ref[1], 2)

    # range over which quantiles are estimated
    P = linspace(0.01, 0.99, rankn)

    #Get date vectors
    datevec_obs = obs[1][Axis{:time}][:]
    datevec_ref = ref[1][Axis{:time}][:]

    # Modify dates (e.g. 29th feb are dropped/lost by default)
    obsvec2, obs_jul, datevec_obs2 = corrjuliandays(obs[1][1,1,:].data, datevec_obs)
    refvec2, ref_jul, datevec_ref2 = corrjuliandays(ref[1][1,1,:].data, datevec_ref)

    # Number of points to sample
    n = round(Int,partition * size(obs[1], 1) * size(obs[1], 2)) # Number of points
    if n == 0
        n=1
    else
        n=n
    end
    # Coordinates of the sampled points
    x = rand(1:size(obs[1],1),n)
    y = rand(1:size(obs[1],2),n)
    # Initialization of the output
    ITP = Array{Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Interpolations.Gridded{typeof(interp)},Tuple{Array{Float64,1}},0},Interpolations.Gridded{typeof(interp)},Interpolations.OnGrid,typeof(extrap)}}(365)
    p = Progress(365, 1)
    # Loop over every julian days
    for ijulian = 1:365
        # Index of ijulian ± window
        idxobs = find_julianday_idx(obs_jul, ijulian, window)
        idxref = find_julianday_idx(ref_jul, ijulian, window)
        # Object containing observation/reference data of the n points on ijulian day
        obsval = fill(NaN, sum(idxobs) * n)
        refval = fill(NaN, sum(idxref) * n)
        for ipoint = 1:n
            iobsvec2, iobs_jul, idatevec_obs2 = corrjuliandays(obs[1][x[ipoint],y[ipoint],:].data, datevec_obs)
            irefvec2, iref_jul, idatevec_ref2 = corrjuliandays(ref[1][x[ipoint],y[ipoint],:].data, datevec_ref)
            obsval[sum(idxobs)*(ipoint-1)+1:sum(idxobs)*ipoint] = iobsvec2[idxobs]
            refval[sum(idxref)*(ipoint-1)+1:sum(idxref)*ipoint] = irefvec2[idxref]
        end

        # Estimate quantiles for obs and ref for ijulian
        obsP = quantile(obsval[.!isnan.(obsval)], P)
        refP = quantile(refval[.!isnan.(refval)], P)
        if lowercase(method) == "additive" # used for temperature
            sf_refP = obsP - refP
        elseif lowercase(method) == "multiplicative" # used for precipitation
            sf_refP = obsP ./ refP
            sf_refP[sf_refP .< 0] = 0.
        end
        # transfert function for ijulian
        itp = interpolate((refP,), sf_refP, Gridded(interp))
        itp = extrapolate(itp, extrap) # add extrapolation
        ITP[ijulian] = itp
        next!(p)
    end
    return ITP
end

"""
    qqmap(fut::ClimGrid, ITP; method::String="Additive")

Quantile-Quantile mapping bias correction with a known transfert function. For each julian day of the year, use the right transfert function to correct model values.

**Options**

**method::String = "Additive" (default) or "Multiplicative"**. Additive is used for most climate variables. Multiplicative is usually bounded variables such as precipitation and humidity.

"""

function qqmap(fut::ClimGrid, ITP; method::String="Additive")
    # Get date vectors
    datevec_fut = fut[1][Axis{:time}][:]
    futvec2, fut_jul, datevec_fut2 = corrjuliandays(fut[1][1,1,:].data, datevec_fut)
    # Prepare output array
    dataout = fill(NaN, (size(fut[1], 1), size(fut[1],2), size(futvec2, 1)))::Array{N, T} where N where T
    # Progress meters
    p = Progress(365, 1)
    # Loop over every points
    # for k = 1:size(fut[1], 2)
    #     for j = 1:size(fut[1], 1)
    #         futvec2, fut_jul, datevec_fut2 = corrjuliandays(fut[1][j,k,:].data, datevec_fut)
            # futvec_corr = similar(futvec2, (size(futvec2)))
            # Loop over every julian day
            for ijulian = 1:365
                idxfut = (fut_jul .== ijulian)
                # Value to correct
                # futval = futvec2[idxfut]
                futval = fut[1][:,:,idxfut].data
                # Transfert function for ijulian
                itp = ITP[ijulian]
                # Correct futval
                if lowercase(method) == "additive" # used for temperature
                    futnew = itp[futval] .+ futval
                elseif lowercase(method) == "multiplicative" # used for precipitation
                    futnew = itp[futval] .* futval
                else
                    error("Wrong method")
                end
                # futvec_corr[idxfut] = futnew
                dataout[:,:,idxfut] = futnew
                next!(p)
            end
            # dataout[j,k,:] = futvec_corr
    #     end
    #     next!(p)
    # end
    lonsymbol = Symbol(fut.dimension_dict["lon"])
    latsymbol = Symbol(fut.dimension_dict["lat"])

    dataout2 = AxisArray(dataout, Axis{lonsymbol}(fut[1][Axis{lonsymbol}][:]), Axis{latsymbol}(fut[1][Axis{latsymbol}][:]),Axis{:time}(datevec_fut2))

    return ClimGrid(dataout2, longrid=fut.longrid, latgrid=fut.latgrid, model=fut.model, project=fut.project, institute=fut.institute, experiment=fut.experiment, run=fut.run, filename=fut.filename, dataunits=fut.dataunits, latunits=fut.latunits, lonunits=fut.lonunits, variable=fut.variable, typeofvar=fut.typeofvar, typeofcal=fut.typeofcal, varattribs=fut.varattribs, globalattribs=fut.varattribs)
end

# function corrjuliandays(obsvec, refvec, futvec, datevec_obs, datevec_ref, datevec_fut)
#
#     # Eliminate February 29th (small price to pay for simplicity and does not affect significantly quantile estimations)
#
#     obs29thfeb = (Dates.month.(datevec_obs) .== Dates.month(Date(2000, 2, 2))) .& (Dates.day.(datevec_obs) .== Dates.day(29))
#     ref29thfeb = (Dates.month.(datevec_ref) .== Dates.month(Date(2000, 2, 2))) .& (Dates.day.(datevec_ref) .== Dates.day(29))
#     fut29thfeb = (Dates.month.(datevec_fut) .== Dates.month(Date(2000, 2, 2))) .& (Dates.day.(datevec_fut) .== Dates.day(29))
#
#     obs_jul = Dates.dayofyear.(datevec_obs)
#     ref_jul = Dates.dayofyear.(datevec_ref)
#     fut_jul = Dates.dayofyear.(datevec_fut)
#
#     # identify leap years
#     leapyears_obs = leapyears(datevec_obs)
#     leapyears_ref = leapyears(datevec_ref)
#     leapyears_fut = leapyears(datevec_fut)
#
#
#     if sum(obs29thfeb) >= 1 & sum(ref29thfeb) == 0 # obs leap year but not models
#
#         for iyear in leapyears_obs
#             k = findfirst(Dates.year.(datevec_obs), iyear) + 59
#             obs_jul[k:k+306] -= 1
#         end
#
#         for iyear in leapyears_ref
#             k = findfirst(Dates.year.(datevec_ref), iyear) + 59
#             ref_jul[k:k+305] -= 1
#         end
#
#         for iyear in leapyears_fut
#             k = findfirst(Dates.year.(datevec_fut), iyear) + 59
#             fut_jul[k:k+305] -= 1
#         end
#
#         datevec_obs2 = datevec_obs[.!obs29thfeb]
#         obsvec2 = obsvec[.!obs29thfeb]
#         obs_jul = obs_jul[.!obs29thfeb]
#
#         refvec2 = refvec
#         datevec_ref2 = datevec_ref
#         futvec2 = futvec
#         datevec_fut2 = datevec_fut
#
#
#         # modify obs_jul to "-=1" for k:k+306 for leap years
#         # modify models ref_jul/fut_jul to "-= 1" for k:k+305 for leap years
#         # delete only obs 29th values
#
#     elseif sum(obs29thfeb) >=1 & sum(ref29thfeb) >=1 # leap years for obs & models
#
#         # modify models obs_jul/ref_jul/fut_jul to "-= 1" for k:k+306 for leap years
#         # delete obs/ref/fut 29th values
#         for iyear in leapyears_obs
#             k = findfirst(Dates.year.(datevec_obs), iyear) + 59
#             obs_jul[k:k+306] -= 1
#         end
#
#         for iyear in leapyears_ref
#             k = findfirst(Dates.year.(datevec_ref), iyear) + 59
#             ref_jul[k:k+306] -= 1
#         end
#
#         for iyear in leapyears_fut
#             k = findfirst(Dates.year.(datevec_fut), iyear) + 59
#             fut_jul[k:k+306] -= 1
#         end
#
#         datevec_obs2 = datevec_obs[.!obs29thfeb]
#         obsvec2 = obsvec[.!obs29thfeb]
#         obs_jul = obs_jul[.!obs29thfeb]
#
#         datevec_ref2 = datevec_ref[.!ref29thfeb]
#         refvec2 = refvec[.!ref29thfeb]
#         ref_jul = ref_jul[.!ref29thfeb]
#
#         datevec_fut2 = datevec_fut[.!fut29thfeb]
#         futvec2 = futvec[.!fut29thfeb]
#         fut_jul = fut_jul[.!fut29thfeb]
#
#     elseif sum(obs29thfeb) == 0 & sum(ref29thfeb) >= 1
#
#         # modify obs_jul to "-=1" for k:k+305 for leap years
#         # modify ref_jul/fut_jul to "-=1" for k:k+306 for leap years
#         # delete ref/fut 29th values
#
#         for iyear in leapyears_obs
#             k = findfirst(Dates.year.(datevec_obs), iyear) + 59
#             obs_jul[k:k+305] -= 1
#         end
#
#         for iyear in leapyears_ref
#             k = findfirst(Dates.year.(datevec_ref), iyear) + 59
#             ref_jul[k:k+306] -= 1
#         end
#
#         for iyear in leapyears_fut
#             k = findfirst(Dates.year.(datevec_fut), iyear) + 59
#             fut_jul[k:k+306] -= 1
#         end
#
#         datevec_obs2 = datevec_obs[.!obs29thfeb]
#         obsvec2 = obsvec[.!obs29thfeb]
#         # obs_jul = obs_jul[.!obs29thfeb]
#
#         datevec_ref2 = datevec_ref[.!ref29thfeb]
#         refvec2 = refvec[.!ref29thfeb]
#         ref_jul = ref_jul[.!ref29thfeb]
#
#         datevec_fut2 = datevec_fut[.!fut29thfeb]
#         futvec2 = futvec[.!fut29thfeb]
#         fut_jul = fut_jul[.!fut29thfeb]
#
#     elseif sum(obs29thfeb) == 0 & sum(ref29thfeb) == 0 # no leap years
#
#         # modify obs_jul/ref_jul/fut_jul to "-=1" for k:k+305 for leap years
#         for iyear in leapyears_obs
#             k = findfirst(Dates.year.(datevec_obs), iyear) + 59
#             obs_jul[k:k+305] -= 1
#         end
#
#         for iyear in leapyears_ref
#             k = findfirst(Dates.year.(datevec_ref), iyear) + 59
#             ref_jul[k:k+305] -= 1
#         end
#
#         for iyear in leapyears_fut
#             k = findfirst(Dates.year.(datevec_fut), iyear) + 59
#             fut_jul[k:k+305] -= 1
#         end
#
#
#     end
#
#     return obsvec2, refvec2, futvec2, obs_jul, ref_jul, fut_jul, datevec_obs2, datevec_ref2, datevec_fut2
#
# end


function corrjuliandays(data_vec, date_vec)

    # Eliminate February 29th (small price to pay for simplicity and does not affect significantly quantile estimations)

    feb29th = (Dates.month.(date_vec) .== Dates.month(Date(2000, 2, 2))) .& (Dates.day.(date_vec) .== Dates.day(29))

    date_jul = Dates.dayofyear.(date_vec)

    # identify leap years
    leap_years = leapyears(date_vec)

    if sum(feb29th) >= 1 # leapyears

        for iyear in leap_years
            k = findfirst(Dates.year.(date_vec), iyear) + 59
            date_jul[k:k+306] -= 1
        end

        date_vec2 = date_vec[.!feb29th]
        data_vec2 = data_vec[.!feb29th]
        date_jul = date_jul[.!feb29th]

    elseif sum(feb29th) == 0 # not a leapyears

        for iyear in leap_years
            k = findfirst(Dates.year.(date_vec), iyear) + 59
            date_jul[k:k+305] -= 1
        end

        date_vec2 = date_vec[.!feb29th]
        data_vec2 = data_vec[.!feb29th]
        # date_jul = date_jul[.!feb29th]

    end

    return data_vec2, date_jul, date_vec2

end

function leapyears(datevec)

    years = unique(Dates.year.(datevec))
    lyrs = years[Dates.isleapyear.(years)]

    return lyrs

end

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
