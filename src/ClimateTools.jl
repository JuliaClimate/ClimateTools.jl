module ClimateTools

# External modules
using NetCDF
using Shapefile
using AxisArrays

# Exported functions
export windnr, leftorright, inpoly, meshgrid, boxcar3, prcp1, frostdays, summerdays, icingdays, tropicalnights, customthresover, customthresunder, annualmax, annualmin, nc2julia, sumleapyear, buildtimevec, inpolyV, shpextract

# Included files
include("functions.jl")
include("indices.jl")
include("extract.jl")


# TYPES
type ClimGrid
  lat::Array{Float64}
  lon::Array{Float64}
  data::AxisArray{Float64}
  timeV::StepRange{Date,Base.Dates.Day}
  model::String
  experiment::String
  run::String
  filename::String
  dataunits::String
  latunits::String
  lonunits::String


  function ClimGrid(lat, lon, data, timeV, model, experiment, run, filename, dataunits, latunits,lonunits)

    # to-do -> add some checks, permutedims if need be
    new(lat, lon, data, timeV, model, experiment, run, filename, dataunits, latunits, lonunits)

  end
end
# TODO Add show method for ClimGrid
# TODO Refactor indices to use AxisArrays
end #module
