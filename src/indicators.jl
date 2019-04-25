"""
vaporpressure(surface_pressure::ClimGrid, specific_humidity::ClimGrid)

Returns the vapor pressure (vp) (Pa) based on the surface pressure (sp) (Pa) and the specific humidity (q).

``vp = \\frac{q * sp}{q+0.622}``

"""
function vaporpressure(specific_humidity::ClimGrid, surface_pressure::ClimGrid)
    @argcheck surface_pressure[9] == "ps"
    @argcheck specific_humidity[9] == "huss"
    @argcheck surface_pressure[2] == "Pa"

    # Calculate vapor pressure
    vp_arraytmp = (specific_humidity.data .* surface_pressure.data) ./ (specific_humidity.data .+ 0.622)
    vp_array = buildarrayinterface(vp_arraytmp, surface_pressure)

    # Build dictionary for the variable vp
    vp_dict = surface_pressure.varattribs
    vp_dict["standard_name"] = "water_vapor_pressure"
    vp_dict["units"] = "Pa"
    vp_dict["history"] = "Water vapor pressure calculated by a function of surface_pressure and specific humidity"

    # Build ClimGrid object
    return ClimGrid(vp_array, longrid=surface_pressure.longrid, latgrid=surface_pressure.latgrid, msk=surface_pressure.msk, grid_mapping=surface_pressure.grid_mapping, dimension_dict=surface_pressure.dimension_dict, timeattrib=surface_pressure.timeattrib, model=surface_pressure.model, frequency=surface_pressure.frequency, experiment=surface_pressure.experiment, run=surface_pressure.run, project=surface_pressure.project, institute=surface_pressure.institute, filename=surface_pressure.filename, dataunits="Pa", latunits=surface_pressure.latunits, lonunits=surface_pressure.lonunits, variable="vp", typeofvar="vp", typeofcal=surface_pressure.typeofcal, varattribs=vp_dict, globalattribs=surface_pressure.globalattribs)
end


"""
vaporpressure(specific_humidity::ClimGrid, sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)

Returns the vapor pressure (vp) (Pa) estimated with the specific humidity (q), the sea level pressure (psl) (Pa), the orography (orog) (m) and the daily mean temperature (tas) (K). An approximation of the surface pressure is first computed by using the sea level pressure, orography and the daily mean temperature (see [`approx_surfacepressure`](@ref)). Then, vapor pressure is calculated by:

``vp = \\frac{q * sp}{q+0.622}``

"""
function vaporpressure(specific_humidity::ClimGrid, sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)

    @argcheck specific_humidity[9] == "huss"
    @argcheck sealevel_pressure[9] == "psl"
    @argcheck orography[9] == "orog"
    @argcheck daily_temperature[9] == "tas"
    @argcheck sealevel_pressure[2] == "Pa"
    @argcheck orography[2] == "m"
    @argcheck in(daily_temperature[2], ["Celsius", "K", "°C"])

    # Convert to Kelvin if necessery
    if daily_temperature[2] == "Celsius" || daily_temperature[2] == "°C"
        daily_temperature = daily_temperature + 273.15
    end

    # Calculate the estimated surface pressure
    surface_pressure = approx_surfacepressure(sealevel_pressure, orography, daily_temperature)

    # Calculate vapor pressure
    vapor_pressure = vaporpressure(specific_humidity, surface_pressure)

    # Return ClimGrid type containing the vapor pressure
    return vapor_pressure
end

"""
approx_surfacepressure(sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)

Returns the approximated surface pressure (*sp*) (Pa) using sea level pressure (*psl*) (Pa), orography (*orog*) (m), and daily mean temperature (*tas*) (K).

``sp = psl * 10^{x}``

where ``x = \\frac{-orog}{18400 * tas / 273.15} ``

"""
function approx_surfacepressure(sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)
    @argcheck sealevel_pressure[9] == "psl"
    @argcheck orography[9] == "orog"
    @argcheck daily_temperature[9] == "tas"

    # Calculate the estimated surface pressure
    exponent = (-1.0 .* orography.data) ./ (18400.0 .* daily_temperature.data ./ 273.15)
    ps_arraytmp = sealevel_pressure.data .* (10.0.^exponent)
    # ps_arraytmp = [ps_arraytmp][1]Pa
    ps_array = buildarrayinterface(ps_arraytmp, sealevel_pressure)

    # Build dictionary for the variable vp
    ps_dict = sealevel_pressure.varattribs
    ps_dict["standard_name"] = "surface_pressure"
    ps_dict["units"] = "Pa"
    ps_dict["history"] = "Surface pressure estimated with the sealevel pressure, the orography and the daily temperature"

    # Build ClimGrid object
    return ClimGrid(ps_array, longrid=sealevel_pressure.longrid, latgrid=sealevel_pressure.latgrid, msk=sealevel_pressure.msk, grid_mapping=sealevel_pressure.grid_mapping, dimension_dict=sealevel_pressure.dimension_dict, timeattrib=sealevel_pressure.timeattrib, model=sealevel_pressure.model, frequency=sealevel_pressure.frequency, experiment=sealevel_pressure.experiment, run=sealevel_pressure.run, project=sealevel_pressure.project, institute=sealevel_pressure.institute, filename=sealevel_pressure.filename, dataunits="Pa", latunits=sealevel_pressure.latunits, lonunits=sealevel_pressure.lonunits, variable="ps", typeofvar="ps", typeofcal=sealevel_pressure.typeofcal, varattribs=ps_dict, globalattribs=sealevel_pressure.globalattribs)
end

"""
wbgt(diurnal_temperature::ClimGrid, vapor_pressure::ClimGrid)

Returns the simplified wet-bulb global temperature (*wbgt*) (Celsius) calculated using the vapor pressure (Pa) of the day and the estimated mean diurnal temperature (Celsius; temperature between 7:00 (7am) and 17:00 (5pm)).

``wbgt = 0.567 * Tday + 0.00393 * vp + 3.94``

"""
function wbgt(mean_temperature::ClimGrid, vapor_pressure::ClimGrid)

    @argcheck istemperature(mean_temperature)#[2], ["Celsius", "K", "°C"])
    @argcheck vapor_pressure[9] == "vp"
    @argcheck istemperature_units(mean_temperature)#[2], ["Celsius", "K", "°C"])
    @argcheck vapor_pressure[2] == "Pa"

    # Convert to Celsius if necessery
    if mean_temperature[2] == "K"
        mean_temperature = mean_temperature - 273.15
    end

    # Calculate the wbgt
    wbgt_arraytmp = (0.567 .* mean_temperature.data) + (0.00393 .* vapor_pressure.data) .+ 3.94
    # wbgt_arraytmp = [wbgt_arraytmp][1]°C
    wbgt_array = buildarrayinterface(wbgt_arraytmp, mean_temperature)

    # Build dictionary for the variable wbgt
    wbgt_dict = mean_temperature.varattribs
    wbgt_dict["standard_name"] = "simplified_wetbulb_globe_temperature"
    wbgt_dict["units"] = "Celsius"
    wbgt_dict["history"] = "Wet-bulb globe temperature estimated with the vapor pressure and the mean temperature"

    # Build ClimGrid object
    return ClimGrid(wbgt_array, longrid=mean_temperature.longrid, latgrid=mean_temperature.latgrid, msk=mean_temperature.msk, grid_mapping=mean_temperature.grid_mapping, dimension_dict=mean_temperature.dimension_dict, timeattrib=mean_temperature.timeattrib, model=mean_temperature.model, frequency=mean_temperature.frequency, experiment=mean_temperature.experiment, run=mean_temperature.run, project=mean_temperature.project, institute=mean_temperature.institute, filename=mean_temperature.filename, dataunits="Celsius", latunits=mean_temperature.latunits, lonunits=mean_temperature.lonunits, variable="wbgt", typeofvar="wbgt", typeofcal=mean_temperature.typeofcal, varattribs=wbgt_dict, globalattribs=mean_temperature.globalattribs)
end

"""
diurnaltemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid, α::Float64)

Returns an estimation of the diurnal temperature (temperature between 7:00 (7am) and 17:00 (5pm)). The estimation is a linear combination of the daily minimum temperature (temperatureminimum) and daily maximum temperature (temperaturemaximum). The value of α has to be estimated seperatly from observations and depends on the location. The daily max and min must be in the same unit and in Celsius or Kelvin The diurnal temperature returned is in the same units as the daily minimum temperature and daily maximum temperature.

``Tdiu = α * Tmin + (1 - α) * Tmax``
"""
function diurnaltemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid, α::Float64)
    @argcheck temperatureminimum[9] == "tasmin"
    @argcheck temperaturemaximum[9] == "tasmax"
    @argcheck istemperature_units(temperatureminimum)#[2], ["Celsius", "K", "°C"])
    @argcheck istemperature_units(temperaturemaximum)#[2], ["Celsius", "K", "°C"])
    @argcheck temperatureminimum[2] == temperaturemaximum[2]

    # Calculate the diurnal temperature
    tdiu = α .*temperatureminimum.data + (1-α) .* temperaturemaximum.data
    # un = get_units(temperatureminimum)
    # tdiu = [tdiu][1]un
    tdiu_array = buildarrayinterface(tdiu, temperatureminimum)

    # Build dictionary for the variable wbgt
    tdiu_dict = temperatureminimum.varattribs
    tdiu_dict["standard_name"] = "diurnal_temperature"
    tdiu_dict["units"] = temperatureminimum[2]
    tdiu_dict["history"] = "Diurnal temperature (between 7:00 and 17:00) estimated using a linear combination of the daily minimum temperature and the daily maximum temperature"

    # Build ClimGrid object
    return ClimGrid(tdiu_array, longrid=temperatureminimum.longrid, latgrid=temperatureminimum.latgrid, msk=temperatureminimum.msk, grid_mapping=temperatureminimum.grid_mapping, dimension_dict=temperatureminimum.dimension_dict, timeattrib=temperatureminimum.timeattrib, model=temperatureminimum.model, frequency=temperatureminimum.frequency, experiment=temperatureminimum.experiment, run=temperatureminimum.run, project=temperatureminimum.project, institute=temperatureminimum.institute, filename=temperatureminimum.filename, dataunits=temperatureminimum.dataunits, latunits=temperatureminimum.latunits, lonunits=temperatureminimum.lonunits, variable="tdiu", typeofvar="tdiu", typeofcal=temperatureminimum.typeofcal, varattribs=tdiu_dict, globalattribs=temperatureminimum.globalattribs)
end

"""
meantemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid)

Returns the daily mean temperature calculated from the maximum and minimum temperature. Daily maximum and minimum temperature must be in the same units. The mean temperature returned is in the same units as the daily minimum temperature and daily maximum temperature.

``Tmean = \\frac{Tmax + Tmin}{2}``
"""
function meantemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid)
    @argcheck temperatureminimum[9] == "tasmin"
    @argcheck temperaturemaximum[9] == "tasmax"
    @argcheck istemperature_units(temperatureminimum)#[2], ["Celsius", "K", "°C"])
    @argcheck istemperature_units(temperaturemaximum)#[2], ["Celsius", "K", "°C"])
    @argcheck temperatureminimum[2] == temperaturemaximum[2]

    # Calculate the diurnal temperature
    # un = unit(temperatureminimum.data[1])
    tmean = (temperatureminimum.data .+ temperaturemaximum.data) ./ 2
    # tmean = [tmean][1]un
    tmean_array = buildarrayinterface(tmean, temperatureminimum)

    # Build dictionary for the variable wbgt
    tmean_dict = temperatureminimum.varattribs
    tmean_dict["standard_name"] = "air_temperature"
    tmean_dict["long_name"] = "mean_daily_temperature"
    tmean_dict["units"] = temperatureminimum[2]
    tmean_dict["history"] = "Mean temperature calculated from the daily minimum and maximum temperature"

    # Build ClimGrid object
    return ClimGrid(tmean_array, longrid=temperatureminimum.longrid, latgrid=temperatureminimum.latgrid, msk=temperatureminimum.msk, grid_mapping=temperatureminimum.grid_mapping, dimension_dict=temperatureminimum.dimension_dict, timeattrib=temperatureminimum.timeattrib, model=temperatureminimum.model, frequency=temperatureminimum.frequency, experiment=temperatureminimum.experiment, run=temperatureminimum.run, project=temperatureminimum.project, institute=temperatureminimum.institute, filename=temperatureminimum.filename, dataunits=temperatureminimum.dataunits, latunits=temperatureminimum.latunits, lonunits=temperatureminimum.lonunits, variable="tmean", typeofvar="tmean", typeofcal=temperatureminimum.typeofcal, varattribs=tmean_dict, globalattribs=temperatureminimum.globalattribs)
end

"""
drought_dc(pr::ClimGrid, tas::Climgrid)

Returns the drought index. The index is defined as a function of daily temperature and daily precipitation. This indice correspond to the DC indice in Wang et al. 2015.

Reference: Wang et al. 2015. Updated source code for calculating fire danger indices in the Canadian Forest Fire Weather Index System. Information Report NOR-X-424. Canadian Forst Service. 36pp.
"""
function drought_dc(pr::ClimGrid, tas::ClimGrid)

    @argcheck ClimateTools.isprecipitation(pr)
    @argcheck ClimateTools.istemperature(tas)
    @argcheck ClimateTools.istemperature_units(tas)
    @argcheck size(pr[1]) == size(tas[1])

    # Convert to Celsius if necessery
    if tas[2] == "K"
        tas = tas - 273.15
    end

    # Check for precip units (needs mm)
    if pr[2] != "mm"
        error("Precipitation needs to be in mm")
    end

    # Time vector
    timevec = get_timevec(pr)

    # Allocated output array
    dataout  = fill(NaN, (size(pr.data, 1), size(pr.data, 2), size(pr.data,3)))

    # Reshape. Allows multi-threading on grid points.
    prin = reshape(pr[1].data, (size(pr[1].data, 1)*size(pr[1].data, 2), size(pr[1].data, 3)))
    tasin = reshape(tas[1].data, (size(tas[1].data, 1)*size(tas[1].data, 2), size(tas[1].data, 3)))
    dataoutin = reshape(dataout, (size(dataout, 1)*size(dataout, 2), size(dataout, 3)))

    # Looping over grid points using multiple-dispatch calls to qqmap
    Threads.@threads for k = 1:size(prin, 1)
    # for k = 1:size(prin, 1)

        prvec = prin[k,:]
        tasvec = tasin[k,:]

        dataoutin[k, :] = drought_dc(prvec, tasvec, timevec)

        # println(k)

    end

    # Build metadata
    dc_dict = pr.varattribs
    dc_dict["standard_name"] = "drought_code"
    dc_dict["units"] = ""
    dc_dict["history"] = "Drought code inspired by Wang et al. 2015. Updated source code for calculating fire danger indices in the Canadian Forest Fire Weather Index System. Information Report NOR-X-424. Canadian Forst Service. 36pp."

    # Build AxisArray
    axisarray = buildarrayinterface(dataout, pr)

    # Build ClimGrid
    return ClimGrid(axisarray, longrid=pr.longrid, latgrid=pr.latgrid, msk=pr.msk, grid_mapping=pr.grid_mapping, dimension_dict=pr.dimension_dict, timeattrib=pr.timeattrib, model=pr.model, frequency=pr.frequency, experiment=pr.experiment, run=pr.run, project=pr.project, institute=pr.institute, filename=pr.filename, dataunits="", latunits=pr.latunits, lonunits=pr.lonunits, variable="dc", typeofvar="dc", typeofcal=pr.typeofcal, varattribs=dc_dict, globalattribs=pr.globalattribs)

end

"""
drought_dc(prvec::Array{N,1}, tasvec::Array{N,1}, timevec)

Returns the drought index. The index is defined as a function of daily temperature and daily precipitation. This indice correspond to the DC indice in Wang et al. 2015.

Reference: Wang et al. 2015. Updated source code for calculating fire danger indices in the Canadian Forest Fire Weather Index System. Information Report NOR-X-424. Canadian Forst Service. 36pp.
"""
function drought_dc(prvec::Array{N,1} where N, tasvec::Array{N,1} where N, timevec)

    dc0 = 15.0 # initial drought code value for each year

    years    = Dates.year.(timevec)
    numYears = unique(years)

    # Output vector
    dc = fill(NaN, size(prvec))

    for i in 1:length(numYears)
        idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])

        # Find index where tas is > 6°C for more than 3 days
        idx_tas = tasvec[idx] .>= 6.0

        nanfrac = sum(idx_tas)#/length(idx_tas)

        if nanfrac >= 3 # need a start and end to the season

            cum_ = cumsum(idx_tas)
            cum_2 = (cum_[1+3:end] - cum_[1:end-3]) .== 3

            if sum(cum_2) >= 2

                idx_start = findfirst(cum_2)
                idx_end = findlast(cum_2)

                # # idx_tas = findall(x -> x >= 6.0, tasvec[idx])
                # idx_start, idx_end = ((findall((diff(idx_tas) .- 1) .== 0) .+ 1)[2],  (findall((diff(idx_tas) .- 1) .== 0) .+ 1)[end])
                # #(findall((diff(idx_tas) .- 1) .== 0) .+ 1)[2]

                # dates_tuple = Dates.yearmonthday.(timevec[idx])

                dc[idx_start] = dc0 # initialize drought code

                for iday = idx_start:idx_end#length(dates_tuple)

                    pr_day = prvec[iday]
                    tas_day = tasvec[iday]

                    # dc₀_temp = dc[iday-1]

                    mth = Dates.month(timevec[iday])

                    dc[iday] = drought_dc(pr_day, tas_day, dc[iday], mth)

                    dc[iday + 1] = dc[iday]

                    # println(iday)


                end
            end
        end

    end

    return dc

end

function drought_dc(pr_day::Real, tas_day::Real, dc0::Real, mth::Int)

    # Months constant
    fl = [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6]
    dc_single = NaN # initialize

    if tas_day < -2.8
        tas_day = -2.8
    end

    pe = (0.36*(tas_day + 2.8) + fl[mth]) / 2

    if pe <= 0.0
        pe = 0.0
    end

    if pr_day > 2.8
        rw = 0.83*pr_day - 1.27
        smi = 800.0*exp(-dc0/400.0)
        if smi == 0.0
            smi = 1.0
        end
        dr = dc0 - 400.0*log(1.0+((3.937*rw/smi)))

        if dr > 0.0
            dc_single = dr + pe
        end

    elseif pr_day <= 2.8
        dc_single = dc0 + pe
    end

    return dc_single


end

function istemperature_units(C::ClimGrid)
    return in(C[2], ["K", "Kelvin", "Celsius", "°C", "C"])
end

function istemperature(C::ClimGrid)
    return in(C[9], ["tasmax", "tas", "tasmin", "tmean", "tmoy"])
end

function isprecipitation(C::ClimGrid)
    return in(C[9], ["pr", "precip", "precipitation", "precipitation_flux"])
end

function isprecipitation_units(C::ClimGrid)
    return in(C[2], ["mm", "kg m-2 s-1", "kg m-2"])
end
