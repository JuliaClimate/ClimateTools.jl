
"""
    vaporpressure(surface_pressure::ClimGrid, specific_humidity::ClimGrid)

Calculation of the vapor pressure (Pa) based on the surface pressure (Pa) and the specific humidity.

Let data[i,j] be the vaporpressure for day i in year j.
"""
function vaporpressure(surface_pressure::ClimGrid, specific_humidity::ClimGrid)
  @argcheck surface_pressure[9] == "ps"
  @argcheck specific_humidity[9] == "huss"

  # Calculate vapor pressure
  vapor_pressure = (specific_humidity * surface_pressure) / (0.622 + specific_humidity)
  vapor_pressure[9] = "vp"

  # Return climGrid type containing the the vapor pressure
  return vapor_pressure
end
