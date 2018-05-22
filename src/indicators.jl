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

Let data[i,j] be the vaporpressure for day i in year j.
"""
function vaporpressure(surface_pressure::ClimGrid, specific_humidity::ClimGrid)
  @argcheck surface_pressure[9] == "ps"
  @argcheck specific_humidity[9] == "huss"
  years    = Dates.year.(surface_pressure.data[Axis{:time}][:])
  numYears = unique(years)
  # Calculate vapor pressure
  vp_arraytmp = (specific_humidity.data .* surface_pressure.data) ./ (specific_humidity.data + 0.622)

  vp_array = buildarrayindicators(vp_arraytmp, surface_pressure)

  vapor_pressure = ClimGrid(vp_array, longrid=surface_pressure.longrid, latgrid=surface_pressure.latgrid, msk=surface_pressure.msk, grid_mapping=surface_pressure.grid_mapping, dimension_dict=surface_pressure.dimension_dict, model=surface_pressure.model, frequency=surface_pressure.frequency, experiment=surface_pressure.experiment, run=surface_pressure.run, project=surface_pressure.project, institute=surface_pressure.institute, filename=surface_pressure.filename, dataunits="Pa", latunits=surface_pressure.latunits, lonunits=surface_pressure.lonunits, variable="vp", typeofvar="vp", typeofcal=surface_pressure.typeofcal, varattribs=surface_pressure.varattribs, globalattribs=surface_pressure.globalattribs)

  # Return climGrid type containing the the vapor pressure
  return vapor_pressure
end
