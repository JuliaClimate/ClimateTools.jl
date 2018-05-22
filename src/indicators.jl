"""
    vaporpressure(surface_pressure::ClimGrid, specific_humidity::ClimGrid)

Calculation of the vapor pressure (vp) (Pa) based on the surface pressure (sp) (Pa) and the specific humidity (q).

``vp = \\frac{q * sp}{q+0.622}``

"""
function vaporpressure(specific_humidity::ClimGrid, surface_pressure::ClimGrid)
  @argcheck surface_pressure[9] == "ps"
  @argcheck specific_humidity[9] == "huss"

  # Calculate vapor pressure
  vp_arraytmp = (specific_humidity.data .* surface_pressure.data) ./ (specific_humidity.data .+ 0.622)
  vp_array = buildarrayinterface(vp_arraytmp, surface_pressure)

  # Build dictionary for the variable vp
  vp_dict = surface_pressure.varattribs
  vp_dict["standard_name"] = "water_vapor_pressure"
  vp_dict["units"] = "Pa"
  vp_dict["units"] = "Water vapor pressure calculated by a function of surface_pressure and specific humidity"

  # Build ClimGrid object
  return ClimGrid(vp_array, longrid=surface_pressure.longrid, latgrid=surface_pressure.latgrid, msk=surface_pressure.msk, grid_mapping=surface_pressure.grid_mapping, dimension_dict=surface_pressure.dimension_dict, model=surface_pressure.model, frequency=surface_pressure.frequency, experiment=surface_pressure.experiment, run=surface_pressure.run, project=surface_pressure.project, institute=surface_pressure.institute, filename=surface_pressure.filename, dataunits="Pa", latunits=surface_pressure.latunits, lonunits=surface_pressure.lonunits, variable="vp", typeofvar="vp", typeofcal=surface_pressure.typeofcal, varattribs=vp_dict, globalattribs=surface_pressure.globalattribs)
end


"""
    vaporpressure(specific_humidity::ClimGrid, sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)

Calculation of the vapor pressure (Pa) estimated with the specific humidity, the sea level pressure (Pa), the orography (m) and the daily temperature (K).

"""
function vaporpressure(specific_humidity::ClimGrid, sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)
  @argcheck specific_humidity[9] == "huss"
  @argcheck sealevel_pressure[9] == "psl"
  @argcheck orography[9] == "orog"
  @argcheck daily_temperature[9] == "tas"

  # Calculate the estimated surface pressure
  surface_pressure = approxim_surfacepressure(sealevel_pressure, orography, daily_temperature)

  # Calculate vapor pressure
  vapor_pressure = vaporpressure(specific_humidity, surface_pressure)

  # Return ClimGrid type containing the vapor pressure
  return vapor_pressure
end

  """
      approxim_surfacepressure(sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)

Calculation of the surface pressure (Pa) approximated from the sea level pressure (Pa), the orography (m), and the daily temperature (K).

"""
function approx_surfacepressure(sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)
  @argcheck sealevel_pressure[9] == "psl"
  @argcheck orography[9] == "orog"
  @argcheck daily_temperature[9] == "tas"

  # Calculate the estimated surface pressure
  exponent = (-1.0 .* orography.data) ./ (18400.0 .* daily_temperature.data ./ 273.15)
  ps_arraytmp = sealevel_pressure.data .* (10.^exponent)
  ps_array = buildarrayinterface(ps_arraytmp, sealevel_pressure)

  # Build dictionary for the variable vp
  ps_dict = sealevel_pressure.varattribs
  ps_dict["standard_name"] = "surface_pressure"
  ps_dict["units"] = "Pa"
  ps_dict["units"] = "Surface pressure estimated with the sealevel pressure, the orography and the daily temperature"

  # Build ClimGrid object
  return ClimGrid(ps_array, longrid=sealevel_pressure.longrid, latgrid=sealevel_pressure.latgrid, msk=sealevel_pressure.msk, grid_mapping=sealevel_pressure.grid_mapping, dimension_dict=sealevel_pressure.dimension_dict, model=sealevel_pressure.model, frequency=sealevel_pressure.frequency, experiment=sealevel_pressure.experiment, run=sealevel_pressure.run, project=sealevel_pressure.project, institute=sealevel_pressure.institute, filename=sealevel_pressure.filename, dataunits="Pa", latunits=sealevel_pressure.latunits, lonunits=sealevel_pressure.lonunits, variable="ps", typeofvar="ps", typeofcal=sealevel_pressure.typeofcal, varattribs=ps_dict, globalattribs=sealevel_pressure.globalattribs)
end
