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
  return ClimGrid(vp_array, longrid=surface_pressure.longrid, latgrid=surface_pressure.latgrid, msk=surface_pressure.msk, grid_mapping=surface_pressure.grid_mapping, dimension_dict=surface_pressure.dimension_dict, model=surface_pressure.model, frequency=surface_pressure.frequency, experiment=surface_pressure.experiment, run=surface_pressure.run, project=surface_pressure.project, institute=surface_pressure.institute, filename=surface_pressure.filename, dataunits="Pa", latunits=surface_pressure.latunits, lonunits=surface_pressure.lonunits, variable="vp", typeofvar="vp", typeofcal=surface_pressure.typeofcal, varattribs=vp_dict, globalattribs=surface_pressure.globalattribs)
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
  @argcheck in(daily_temperature[2], ["Celsius", "K"])

  # Convert to Kelvin if necessery
  if daily_temperature[2] == "Celsius"
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
  ps_array = buildarrayinterface(ps_arraytmp, sealevel_pressure)

  # Build dictionary for the variable vp
  ps_dict = sealevel_pressure.varattribs
  ps_dict["standard_name"] = "surface_pressure"
  ps_dict["units"] = "Pa"
  ps_dict["history"] = "Surface pressure estimated with the sealevel pressure, the orography and the daily temperature"

  # Build ClimGrid object
  return ClimGrid(ps_array, longrid=sealevel_pressure.longrid, latgrid=sealevel_pressure.latgrid, msk=sealevel_pressure.msk, grid_mapping=sealevel_pressure.grid_mapping, dimension_dict=sealevel_pressure.dimension_dict, model=sealevel_pressure.model, frequency=sealevel_pressure.frequency, experiment=sealevel_pressure.experiment, run=sealevel_pressure.run, project=sealevel_pressure.project, institute=sealevel_pressure.institute, filename=sealevel_pressure.filename, dataunits="Pa", latunits=sealevel_pressure.latunits, lonunits=sealevel_pressure.lonunits, variable="ps", typeofvar="ps", typeofcal=sealevel_pressure.typeofcal, varattribs=ps_dict, globalattribs=sealevel_pressure.globalattribs)
end

"""
    wbgt(diurnal_temperature::ClimGrid, vapor_pressure::ClimGrid)

Returns the simplified wet-bulb global temperature (*wbgt*) (Celsius) calculated using the vapor pressure (Pa) of the day and the estimated mean diurnal temperature (Celsius; temperature between 7:00 (7am) and 17:00 (5pm)).

``wbgt = 0.567 * Tday + 0.00393 * vp + 3.94``

"""
function wbgt(mean_temperature::ClimGrid, vapor_pressure::ClimGrid)
  @argcheck in(mean_temperature[9], ["tmean", "tasmax", "tasmin"])
  @argcheck vapor_pressure[9] == "vp"
  @argcheck in(mean_temperature[2], ["Celsius", "K"])
  @argcheck vapor_pressure[2] == "Pa"

  # Convert to Celsius if necessery
  if mean_temperature[2] == "K"
    mean_temperature -= 273.15
  end

  # Calculate the wbgt
  wbgt_arraytmp = (0.567 .* mean_temperature.data) + (0.00393 .* vapor_pressure.data) .+ 3.94
  wbgt_array = buildarrayinterface(wbgt_arraytmp, mean_temperature)

  # Build dictionary for the variable wbgt
  wbgt_dict = mean_temperature.varattribs
  wbgt_dict["standard_name"] = "simplified_wetbulb_globe_temperature"
  wbgt_dict["units"] = "Celsius"
  wbgt_dict["history"] = "Wet-bulb globe temperature estimated with the vapor pressure and the mean temperature"

  # Build ClimGrid object
  return ClimGrid(wbgt_array, longrid=mean_temperature.longrid, latgrid=mean_temperature.latgrid, msk=mean_temperature.msk, grid_mapping=mean_temperature.grid_mapping, dimension_dict=mean_temperature.dimension_dict, model=mean_temperature.model, frequency=mean_temperature.frequency, experiment=mean_temperature.experiment, run=mean_temperature.run, project=mean_temperature.project, institute=mean_temperature.institute, filename=mean_temperature.filename, dataunits="Celsius", latunits=mean_temperature.latunits, lonunits=mean_temperature.lonunits, variable="wbgt", typeofvar="wbgt", typeofcal=mean_temperature.typeofcal, varattribs=wbgt_dict, globalattribs=mean_temperature.globalattribs)
end

"""
    diurnaltemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid, α::Float64)

Returns an estimation of the diurnal temperature (temperature between 7:00 (7am) and 17:00 (5pm)). The estimation is a linear combination of the daily minimum temperature (temperatureminimum) and daily maximum temperature (temperaturemaximum). The value of α has to be estimated seperatly from observations and depends on the location. The daily max and min must be in the same unit and in Celsius or Kelvin The diurnal temperature returned is in the same units as the daily minimum temperature and daily maximum temperature.

``Tdiu = α * Tmin + (1 - α) * Tmax``
"""
function diurnaltemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid, α::Float64)
  @argcheck temperatureminimum[9] == "tasmin"
  @argcheck temperaturemaximum[9] == "tasmax"
  @argcheck in(temperatureminimum[2], ["Celsius", "K"])
  @argcheck in(temperaturemaximum[2], ["Celsius", "K"])
  @argcheck temperatureminimum[2] == temperaturemaximum[2]

  # Calculate the diurnal temperature
  tdiu = α .*temperatureminimum.data + (1-α) .* temperaturemaximum.data
  tdiu_array = buildarrayinterface(tdiu, temperatureminimum)

  # Build dictionary for the variable wbgt
  tdiu_dict = temperatureminimum.varattribs
  tdiu_dict["standard_name"] = "diurnal_temperature"
  tdiu_dict["units"] = temperatureminimum[2]
  tdiu_dict["history"] = "Diurnal temperature (between 7:00 and 17:00) estimated using a linear combination of the daily minimum temperature and the daily maximum temperature"

  # Build ClimGrid object
  return ClimGrid(tdiu_array, longrid=temperatureminimum.longrid, latgrid=temperatureminimum.latgrid, msk=temperatureminimum.msk, grid_mapping=temperatureminimum.grid_mapping, dimension_dict=temperatureminimum.dimension_dict, model=temperatureminimum.model, frequency=temperatureminimum.frequency, experiment=temperatureminimum.experiment, run=temperatureminimum.run, project=temperatureminimum.project, institute=temperatureminimum.institute, filename=temperatureminimum.filename, dataunits=temperatureminimum.dataunits, latunits=temperatureminimum.latunits, lonunits=temperatureminimum.lonunits, variable="tdiu", typeofvar="tdiu", typeofcal=temperatureminimum.typeofcal, varattribs=tdiu_dict, globalattribs=temperatureminimum.globalattribs)
end

"""
    meantemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid)

Returns the daily mean temperature calculated from the maximum and minimum temperature. Daily maximum and minimum temperature must be in the same units. The mean temperature returned is in the same units as the daily minimum temperature and daily maximum temperature.

``Tmean = \\frac{Tmax + Tmin}{2}``
"""
function meantemperature(temperatureminimum::ClimGrid, temperaturemaximum::ClimGrid)
  @argcheck temperatureminimum[9] == "tasmin"
  @argcheck temperaturemaximum[9] == "tasmax"
  @argcheck in(temperatureminimum[2], ["Celsius", "K"])
  @argcheck in(temperaturemaximum[2], ["Celsius", "K"])
  @argcheck temperatureminimum[2] == temperaturemaximum[2]

  # Calculate the diurnal temperature
  tmean = (temperatureminimum.data .+ temperaturemaximum.data) ./ 2
  tmean_array = buildarrayinterface(tmean, temperatureminimum)

  # Build dictionary for the variable wbgt
  tmean_dict = temperatureminimum.varattribs
  tmean_dict["standard_name"] = "air_temperature"
  tmean_dict["long_name"] = "mean_daily_temperature"
  tmean_dict["units"] = temperatureminimum[2]
  tmean_dict["history"] = "Mean temperature calculated from the daily minimum and maximum temperature"

  # Build ClimGrid object
  return ClimGrid(tmean_array, longrid=temperatureminimum.longrid, latgrid=temperatureminimum.latgrid, msk=temperatureminimum.msk, grid_mapping=temperatureminimum.grid_mapping, dimension_dict=temperatureminimum.dimension_dict, model=temperatureminimum.model, frequency=temperatureminimum.frequency, experiment=temperatureminimum.experiment, run=temperatureminimum.run, project=temperatureminimum.project, institute=temperatureminimum.institute, filename=temperatureminimum.filename, dataunits=temperatureminimum.dataunits, latunits=temperatureminimum.latunits, lonunits=temperatureminimum.lonunits, variable="tmean", typeofvar="tmean", typeofcal=temperatureminimum.typeofcal, varattribs=tmean_dict, globalattribs=temperatureminimum.globalattribs)
end
