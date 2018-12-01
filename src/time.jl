"""
    function temporalsubset(C::ClimGrid, startdate::Date, enddate::Date)

Returns the temporal subset of ClimGrid C. The temporal subset is defined by a start and end date.

"""
function temporalsubset(C::ClimGrid, datebeg::Tuple, dateend::Tuple)

    T = typeof(get_timevec(C)[1])
    timeV = get_timevec(C)
    idxtimebeg, idxtimeend = ClimateTools.timeindex(timeV, datebeg, dateend, T)

    startdate = ClimateTools.buildtimetype(datebeg, T)
    enddate = ClimateTools.buildtimetype(dateend, T)

    # some checkups
    @argcheck startdate <= enddate
    # @argcheck startdate >= C[1][Axis{:time}][1]
    # @argcheck enddate <= C[1][Axis{:time}][end]

    dataOut = C[1][Axis{:time}(idxtimebeg:idxtimeend)]

    # The following control ensure that a 1-timestep temporal subset returns a 3D Array with time information on the timestep. i.e. startdate == enddate
    if ndims(dataOut) == 2
        timeV = startdate
        latsymbol = Symbol(C.dimension_dict["lat"])
        lonsymbol = Symbol(C.dimension_dict["lon"])
        data2 = fill(NaN, (size(dataOut,1), size(dataOut, 2), 1))
        data2[:,:,1] = dataOut
        dataOut = AxisArray(data2, Axis{lonsymbol}(C[1][Axis{lonsymbol}][:]), Axis{latsymbol}(C[1][Axis{latsymbol}][:]), Axis{:time}(C[1][Axis{:time}][:]))

    end

    return ClimGrid(dataOut, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)

end

"""
    resample(C::ClimGrid, startmonth::Int64, endmonth::Int64)

Return a resampled subset of ClimGrid C based on months.
"""
function resample(C::ClimGrid, startmonth::Int64, endmonth::Int64)

    @argcheck startmonth >= minimum(Dates.month.(C[1][Axis{:time}][:]))
    @argcheck startmonth <= maximum(Dates.month.(C[1][Axis{:time}][:]))

    if startmonth <= endmonth
        # Each matrix [:,:,i] represent data for a day
        datain = C.data.data
        # Date vector
        datevecin = C[1][Axis{:time}][:]
        # Where are the data between startmonth and endmonth
        index = (Dates.month.(datevecin) .<= endmonth) .&  (Dates.month.(datevecin) .>= startmonth)
        # Keep only data between startmonth and endmonth
        dataout = datain[:,:,index]
        datevecout = datevecin[index]
        # Create the ClimGrid output
        lonsymbol = Symbol(C.dimension_dict["lon"])
        latsymbol = Symbol(C.dimension_dict["lat"])
        axisout = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}][:]), Axis{latsymbol}(C[1][Axis{latsymbol}][:]),    Axis{:time}(datevecout))

    elseif endmonth < startmonth
        # Each matrix [:,:,i] represent data for a day
        datain = C.data.data
        # Date vector
        datevecin = C[1][Axis{:time}][:]
        # Where are the data between startmonth and endmonth
        index = (Dates.month.(datevecin) .<= endmonth) .|  (Dates.month.(datevecin) .>= startmonth)
        # Keep only data between startmonth and endmonth
        dataout = datain[:,:,index]
        datevecout = datevecin[index]
        # Create the ClimGrid output
        lonsymbol = Symbol(C.dimension_dict["lon"])
        latsymbol = Symbol(C.dimension_dict["lat"])
        axisout = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}][:]), Axis{latsymbol}(C[1][Axis{latsymbol}][:]), Axis{:time}(datevecout))
    end

    return ClimGrid(axisout, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping,dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project,institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable,typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    resample(C::ClimGrid, season::String)

Return a resampled subset of ClimGrid C for a given season. Season options are: "DJF" (December-February), "MAM" (March-May), "JJA" (June-August), "SON" (September-November)
"""
function resample(C::ClimGrid, season::String)
    if lowercase(season) == "djf"
      D = resample(C, 12, 2)
    elseif lowercase(season) == "mam"
      D = resample(C, 3, 5)
    elseif lowercase(season) == "jja"
      D = resample(C, 6, 8)
    elseif lowercase(season) == "son"
      D = resample(C, 9, 11)
    else
      error("Wrong season name. Options are DJF (December-February; Winter), MAM (March-May; Spring), JJA (June-August; Summer) and SON (September-November; Fall)")
    end

    return D
end

"""
    periodmean(C::ClimGrid; startdate::Tuple, enddate::Tuple)

Mean of array data over a given period.
"""
function periodmean(C::ClimGrid; start_date::Tuple=(Inf, ), end_date::Tuple=(Inf,))

    if start_date == (Inf, )
        timevec = get_timevec(C)
        # Get time resolution
        rez = C.frequency
        if rez == "year"
            start_date = (Dates.year(timevec[1]), )
            end_date = (Dates.year(timevec[end]), )
        else
            start_date = (Dates.year(timevec[1]), Dates.month(timevec[1]), Dates.day(timevec[1]), Dates.hour(timevec[1]), Dates.minute(timevec[1]), Dates.second(timevec[1]))
            end_date = (Dates.year(timevec[end]), Dates.month(timevec[end]), Dates.day(timevec[end]), Dates.hour(timevec[end]), Dates.minute(timevec[end]), Dates.second(timevec[end]))
        end
    end
    Csubset = temporalsubset(C, start_date, end_date)
    datain   = Csubset.data.data

    # Mean and squeeze
    dataout = fill(NaN, size(datain, 1), size(datain, 2))
    dataout_rshp = reshape(dataout, (size(dataout, 1)*size(dataout, 2)))
    datain_rshp = reshape(datain, (size(datain, 1)*size(datain, 2), size(datain, 3)))

    Threads.@threads for k = 1:size(datain_rshp, 1)
        datatmp = datain_rshp[k, :]
        dataout_rshp[k] = Statistics.mean(datatmp[.!isnan.(datatmp)])
    end

    # Build output AxisArray
    FD = buildarray_climato(C, dataout)

    # Return climGrid type containing the indice
    return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable="periodmean", typeofvar=C.typeofvar, typeofcal="climatology", varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    timeindex(timeVec, start_date, end_date, T)

Return the index of time vector specified by start_date and end_date. T is the DateTime type (see NCDatasets.jl documentation).
"""
function timeindex(timeV, datebeg::Tuple, dateend::Tuple, T)

    # Start Date
    if !isinf(datebeg[1])
        # Build DateTime type
        start_date = ClimateTools.buildtimetype(datebeg, T)
        # @argcheck start_date >= timeV[1]
        idxtimebeg = findfirst(timeV .>= start_date)[1]
    else
        idxtimebeg = 1
    end
    # End date
    if !isinf(dateend[1])
        # Build DateTime type
        end_date = ClimateTools.buildtimetype(dateend, T)
        # @argcheck end_date <= timeV[end]
        idxtimeend = findlast(timeV .<= end_date)[1]
    else
        idxtimeend = length(timeV)
    end

    if !isinf(datebeg[1]) && !isinf(dateend[1])
        @argcheck start_date <= end_date
    end
    return idxtimebeg, idxtimeend
end

# function timeindex(C::ClimGrid, datebeg::Tuple, dateend::Tuple)
#
#     T = typeof(get_timevec(C)[1])
#
#     # Start Date
#     if !isinf(datebeg[1])
#
#         # Build DateTime type
#         start_date = ClimateTools.buildtimetype(datebeg, T)
#
#         # Check
#         @argcheck start_date >= timeV[1]
#
#         idxtimebeg = findfirst(timeV .>= start_date)[1]
#     else
#         idxtimebeg = 1
#     end
#
#     # End date
#     if !isinf(dateend[1])
#
#         # Build DateTime type
#         end_date = ClimateTools.buildtimetype(dateend, T)
#
#         @argcheck end_date <= timeV[end]
#
#         idxtimeend = findlast(timeV .<= end_date)[1]
#     else
#         idxtimeend = length(timeV)
#     end
#
#     if !isinf(datebeg[1]) && !isinf(dateend[1])
#         @argcheck start_date <= end_date
#     end
#
#     return idxtimebeg, idxtimeend
# end

"""
    buildtimetype(datetuple, f)

Returns the adequate DateTime for temporal subsetting using DateType *f*
"""
function buildtimetype(date_tuple, f)

    if length(date_tuple) == 1
        dateout = f(date_tuple[1], 01, 01)
    elseif length(date_tuple) == 2
        dateout = f(date_tuple[1], date_tuple[2], 01)
    elseif length(date_tuple) == 3
        dateout = f(date_tuple[1], date_tuple[2], date_tuple[3])
    elseif length(date_tuple) == 4
        dateout = f(date_tuple[1], date_tuple[2], date_tuple[3], date_tuple[4], 00, 00)
    elseif length(date_tuple) == 5
        dateout = f(date_tuple[1], date_tuple[2], date_tuple[3], date_tuple[4], date_tuple[5], 00)
    elseif length(date_tuple) == 6
        dateout = f(date_tuple[1], date_tuple[2], date_tuple[3], date_tuple[4], date_tuple[5], date_tuple[6])
    end

    return dateout
end

"""
    daymean(C::ClimGrid)

Returns the daily average given a sub-daily ClimGrid.
"""
function daymean(C::ClimGrid)

    datain = C[1].data

    timevec   = get_timevec(C)
    years     = Dates.year.(timevec)
    numYears  = unique(years)
    months    = Dates.month.(timevec)
    numMonths = unique(months)
    days    = Dates.day.(timevec)
    numDays = unique(days)

    T = typeof(timevec[1])

    # dayfactor = ClimateTools.daymean_factor(C.frequency)
    nbdays = length(unique(yearmonthday.(timevec)))
    dataout = zeros(typeof(datain[1]), (size(C[1], 1), size(C[1], 2), nbdays))
    newtime = Array{T}(undef, nbdays)

    # loop over year-month-days
    z = 1
    for iyear in 1:length(numYears)
        for imonth in 1:length(numMonths)
            for iday in 1:daysinmonth(T(numYears[iyear],numMonths[imonth]))

                datefind = T(numYears[iyear],numMonths[imonth],numDays[iday])

                idx = findall(x -> Dates.year(x) == Dates.year(datefind) && Dates.month(x) == Dates.month(datefind) && Dates.day(x) == Dates.day(datefind), timevec)

                dataout[:, :, z] = Statistics.mean(datain[:, :, idx], dims=3)

                newtime[z] = datefind
                z += 1
            end
        end
    end

    # Build output AxisArray
    FD = buildarray_resample(C, dataout, newtime)

    return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="day", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)

end

"""
    daysum(C::ClimGrid)

Returns the daily sum given a sub-daily ClimGrid C.
"""
function daysum(C::ClimGrid)

    datain = C[1].data

    timevec   = get_timevec(C)
    years     = Dates.year.(timevec)
    numYears  = unique(years)
    months    = Dates.month.(timevec)
    numMonths = unique(months)
    days    = Dates.day.(timevec)
    numDays = unique(days)

    T = typeof(timevec[1])

    # dayfactor = ClimateTools.daymean_factor(C.frequency)
    nbdays = length(unique(yearmonthday.(timevec)))
    dataout = zeros(typeof(datain[1]), (size(C[1], 1), size(C[1], 2), nbdays))
    newtime = Array{T}(undef, nbdays)

    # loop over year-month-days
    z = 1
    for iyear in 1:length(numYears)
        for imonth in 1:length(numMonths)
            datemonth = T(numYears[iyear],numMonths[imonth])
            idx_day = findall(x -> Dates.year(x) == Dates.year(datemonth) && Dates.month(x) == Dates.month(datemonth), timevec)
            
            for iday in 1:daysinmonth(T(numYears[iyear],numMonths[imonth]))

                datefind = T(numYears[iyear],numMonths[imonth],numDays[iday])

                idx = findall(x -> Dates.year(x) == Dates.year(datefind) && Dates.month(x) == Dates.month(datefind) && Dates.day(x) == Dates.day(datefind), timevec)

                dataout[:, :, z] = Statistics.sum(datain[:, :, idx], dims=3)

                newtime[z] = datefind
                z += 1
            end
        end
    end

    # Build output AxisArray
    FD = buildarray_resample(C, dataout, newtime)

    return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency="day", experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)

end

# """
#     daysinmonth(D::DateTimeNoLeap)
#
# Workaround to work with non-standard calendar DateTimeNoLeap.
# """
# function daysinmonth(D::DateTimeNoLeap)
#     DAYSINMONTH2 = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#     return DAYSINMONTH2[Dates.month(D)]
# end
#
# """
#     daysinmonth(D::DateTimeStandard)
#
# Workaround to work with non-standard calendar DateTimeStandard.
# """
# function daysinmonth(D::DateTimeStandard)
#     DAYSINMONTH2 = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#     return DAYSINMONTH2[Dates.month(D)] + (Dates.month(D) == 2 && isleapyear(Dates.year(D)))
# end
#
# """
#     daysinmonth(D::DateTimeProlepticGregorian)
#
# Workaround to work with non-standard calendar DateTimeProlepticGregorian.
# """
# function daysinmonth(D::DateTimeProlepticGregorian)
#     DAYSINMONTH2 = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#     return DAYSINMONTH2[m] + (m == 2 && isleapyear(y))
# end
#
# """
#     daysinmonth(D::DateTimeAllLeap)
#
# Workaround to work with non-standard calendar DateTimeAllLeap.
# """
# function daysinmonth(D::DateTimeAllLeap)
#     DAYSINMONTH2 = (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
#     return DAYSINMONTH2[Dates.month(D)]
# end
#
# """
#     daysinmonth(D::DateTime360Day)
#
# Workaround to work with non-standard calendar DateTime360Day.
# """
# daysinmonth(D::DateTime360Day) = 30

"""
    buildtoken(date::DateTime, C::ClimGrid)

Returns the right date token based on time vector contained in ClimGrid C.
"""
function buildtoken(date, C::ClimGrid)
    timevec = get_timevec(C)
    typetoken = typeof(timevec[1])

    return typetoken

end



"""
    timeresolution(timevec::Array{N,1} where N)

Return the time resolution of the vector timevec.

"""
function timeresolution(timevec::Array{N,1} where N)

    # timevec = (NetCDF.ncread(str, "time"))
    if length(timevec) > 1
        timediff = diff(timevec)[2]
        if timediff == 1. || timediff == 1
            return "24h"
        elseif round(timediff, digits=5) == round(12/24, digits=5)
            return "12h"
        elseif round(timediff, digits=5) == round(6/24, digits=5)
            return "6h"
        elseif round(timediff, digits=5) == round(3/24, digits=5)
            return "3h"
        elseif round(timediff, digits=5) == round(1/24, digits=5)
            return "1h"
        elseif round(timediff, digits=5) == 365.0 || round(timediff, digits=5) == 366.0 || round(timediff, digits=5) == 360.0
            return "Yearly"
        end
    else
        return "N/A"
    end
end

"""
    function pr_timefactor(rez::String)

Return the time factor that should be applied to precipitation to get accumulation for resolution "rez"

"""
function pr_timefactor(rez::String)

    if rez == "24h"
        return 86400.0
    elseif rez == "12h"
        return 43200.0
    elseif rez == "6h"
        return 21600.0
    elseif rez == "3h"
        return 10800.0
    elseif rez == "1h"
        return 3600.0
    elseif rez == "N/A"
        return 1.0
    end

end

"""
    function daymean_factor(rez::String)

Return the time factor that should be applied to precipitation to get accumulation for resolution "rez"

"""
function daymean_factor(rez::String)

    if rez == "24h" || rez == "day" || rez == "daily"
        return 1
    elseif rez == "12h"
        return 2
    elseif rez == "6h"
        return 4
    elseif rez == "3h"
        return 8
    elseif rez == "1h"
        return 24
    elseif rez == "N/A"
        return 1
    end

end

"""
    correctdate(C::ClimGrid)

Removes february 29th. Needed for bias correction.
"""
function correctdate(C::ClimGrid)
    date_vec = get_timevec(C) # [1][Axis{:time}][:]
    f = typeof(date_vec[1])
    feb29th = (Dates.month.(date_vec) .== Dates.month(f(2000, 2, 2))) .& (Dates.day.(date_vec) .== Dates.day(29))
    dataout = C[1][:, :, .!feb29th]
    return ClimGrid(dataout; longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    leapyears(datevec)

Returns the leap years contained into datevec.
"""
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

"""
    corrjuliandays(data_vec, date_vec)

Removes 29th february and correct the associated julian days.
"""
function corrjuliandays(data_vec, date_vec)

    # Eliminate February 29th (small price to pay for simplicity and does not affect significantly quantile estimations)

    feb29th = (Dates.month.(date_vec) .== Dates.month(Date(2000, 2, 2))) .& (Dates.day.(date_vec) .== Dates.day(29))

    date_jul = Dates.dayofyear.(date_vec)

    # identify leap years
    leap_years = leapyears(date_vec)

    if sum(feb29th) >= 1 # leapyears

        for iyear in leap_years

            days = date_jul[Dates.year.(date_vec) .== iyear] # days for iyear

            if days[1] >=60 # if the year starts after Feb 29th
                k1 = something(findfirst(isequal(iyear), Dates.year.(date_vec)), 0)
                # k1 = findfirst(Dates.year.(date_vec), iyear) # k1 is the first day
            else
                k1 = something(findfirst(isequal(iyear), Dates.year.(date_vec)), 0) + 60 -days[1]
                # k1 = findfirst(Dates.year.(date_vec), iyear) + 60 - days[1] # else k1 (60-first_julian_day) of the year
            end
            k2 = something(findlast(isequal(iyear), Dates.year.(date_vec)), 0)
            # k2 = findlast(Dates.year.(date_vec), iyear) #+ length(days) - 1 #the end of the year is idx of the first day + number of days in the year - 1
            # k = findfirst(Dates.year.(date_vec), iyear) + 59
            date_jul[k1:k2] .-= 1
        end

        date_vec2 = date_vec[.!feb29th]
        data_vec2 = data_vec[.!feb29th]
        date_jul2 = date_jul[.!feb29th]

    elseif sum(feb29th) == 0 # not a leapyears

        for iyear in leap_years
            days = date_jul[Dates.year.(date_vec) .== iyear] # days for iyear
            if days[1] >= 60 # if the year starts after Feb 29th
                k1 = something(findfirst(isequal(iyear), Dates.year.(date_vec)), 0)
                # k1 = findfirst(Dates.year.(date_vec), iyear) # k1 is the first day
            else
                k1 = something(findfirst(isequal(iyear), Dates.year.(date_vec)), 0) + 60 - days[1]
                # k1 = findfirst(Dates.year.(date_vec), iyear) + 60 - days[1] # else k1 (60-first_julian_day) of the year
            end
            k2 = something(findlast(isequal(iyear), Dates.year.(date_vec)), 0)
            # k2 = findlast(Dates.year.(date_vec), iyear) #+ length(days) - 1 #the end of the year is idx of the first day + number of days in the year - 1
            # k = findfirst(Dates.year.(date_vec), iyear) + 59
            date_jul[k1:k2] .-= 1
        end

        date_vec2 = date_vec[.!feb29th]
        data_vec2 = data_vec[.!feb29th]
        date_jul2 = date_jul[.!feb29th]

    end

    return data_vec2, date_jul2, date_vec2

end
