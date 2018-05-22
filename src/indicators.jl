# Same function as buildarrayinterface
function buildarrayindicators(axisArraytmp, A)
    latsymbol = Symbol(A.dimension_dict["lat"])
    lonsymbol = Symbol(A.dimension_dict["lon"])
    axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}][:]), Axis{latsymbol}(A[1][Axis{latsymbol}][:]), Axis{:time}(A[1][Axis{:time}][:]))
    return axisArray
end

"""
    vaporpressure(surface_pressure::ClimGrid, specific_humidity::ClimGrid)

Calculation of the vapor pressure (Pa) based on the surface pressure (Pa) and the specific humidity.

Let data[i,j] be the vapor pressure for day i in year j.
"""
function vaporpressure(specific_humidity::ClimGrid, surface_pressure::ClimGrid)
  @argcheck surface_pressure[9] == "ps"
  @argcheck specific_humidity[9] == "huss"

  # Calculate vapor pressure
  vp_arraytmp = (specific_humidity.data .* surface_pressure.data) ./ (specific_humidity.data + 0.622)
  vp_array = buildarrayindicators(vp_arraytmp, surface_pressure)

  # Build ClimGrid object
  vapor_pressure = ClimGrid(vp_array, longrid=surface_pressure.longrid, latgrid=surface_pressure.latgrid, msk=surface_pressure.msk, grid_mapping=surface_pressure.grid_mapping, dimension_dict=surface_pressure.dimension_dict, model=surface_pressure.model, frequency=surface_pressure.frequency, experiment=surface_pressure.experiment, run=surface_pressure.run, project=surface_pressure.project, institute=surface_pressure.institute, filename=surface_pressure.filename, dataunits="Pa", latunits=surface_pressure.latunits, lonunits=surface_pressure.lonunits, variable="vp", typeofvar="vp", typeofcal=surface_pressure.typeofcal, varattribs=surface_pressure.varattribs, globalattribs=surface_pressure.globalattribs)

  # Return ClimGrid type containing the vapor pressure
  return vapor_pressure
end


"""
    vaporpressure(specific_humidity::ClimGrid, sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)

Calculation of the vapor pressure (Pa) estimated with the specific humidity, the sea level pressure (Pa), the orography (m) and the daily temperature (K).

Let data[i,j] be the vapor pressure for day i in year j.
"""
function vaporpressure(specific_humidity::ClimGrid, sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)
  @argcheck specific_humidity[9] == "huss"
  @argcheck sealevel_pressure[9] == "psl"
  @argcheck orography[9] == "oro"
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

Let data[i,j] be the surface pressure for day i in year j
"""
function approxim_surfacepressure(sealevel_pressure::ClimGrid, orography::ClimGrid, daily_temperature::ClimGrid)
  @argcheck sealevel_pressure[9] == "psl"
  @argcheck orography[9] == "oro"
  @argcheck daily_temperature[9] == "tas"

  # Calculate the estimated surface pressure
  exponent = (-1*orography) ./ (18400*daily_temperature/273.15)
  ps_arraytmp = sealevel_pressure .* (10.^exponent)
  ps_array = buildarrayindicators(ps_arraytmp, sealevel_pressure)

  # Build ClimGrid object
  surface_pressure = ClimGrid(ps_array, longrid=sealevel_pressure.longrid, latgrid=sealevel_pressure.latgrid, msk=sealevel_pressure.msk, grid_mapping=sealevel_pressure.grid_mapping, dimension_dict=sealevel_pressure.dimension_dict, model=sealevel_pressure.model, frequency=sealevel_pressure.frequency, experiment=sealevel_pressure.experiment, run=sealevel_pressure.run, project=sealevel_pressure.project, institute=sealevel_pressure.institute, filename=sealevel_pressure.filename, dataunits="Pa", latunits=sealevel_pressure.latunits,