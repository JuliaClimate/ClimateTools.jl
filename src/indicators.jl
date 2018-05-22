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
  vp_array = buildarrayinterface(vp_arraytmp, surface_pressure)

  # Build dictionary for the variable vp
  vp_dict = surface_pressure.varattribs
  vp_dict["standard_name"] = "water_vapor_pressure"
  vp_dict["units"] = "Pa"
  vp_dict["units"] = "Water vapor pressure calculated by a function of surface_pressure and specific humidity"

  # Build ClimGrid object
  return ClimGrid(vp_array, longrid=surface_pressure.longrid, latgrid=surface_pressure.latgrid, msk=surface_pressure.msk, grid_mapping=surface_pressure.grid_mapping, dimension_dict=surface_pressure.dimension_dict, model=surface_pressure.model, frequency=surface_pressure.frequency, experiment=surface_pressure.experiment, run=surface_pressure.run, project=surface_pressure.project, institute=surface_pressure.institute, filename=surface_pressure.filename, dataunits="Pa", latunits=surface_pressure.latunits, lonunits=surface_pressure.lonunits, variable="vp", typeofvar="vp", typeofcal=surface_pressure.typeofcal, varattribs=vp_dict, globalattribs=surface_pressure.globalattribs)
end
