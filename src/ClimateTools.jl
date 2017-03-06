module ClimateTools

# External modules
using NetCDF
using Shapefile
using AxisArrays
import Base.vcat

# Exported functions
export windnr, leftorright, inpoly, meshgrid, boxcar3, prcp1, frostdays, summerdays, icingdays, tropicalnights, customthresover, customthresunder, annualmax, annualmin, nc2julia, sumleapyear, buildtimevec, inpolyV, shpextract


# TYPES
type ClimGrid
  data::AxisArray
  model::String
  experiment::String
  run::String
  filename::String
  dataunits::String
  latunits::String
  lonunits::String

  function ClimGrid(data, model = "", experiment = "", run = "", filename = "", dataunits = "", latunits = "",lonunits = "")

    # to-do -> add some checks, permutedims if need be
    new(data, model, experiment, run, filename, dataunits, latunits, lonunits)

  end
end

# Included files
include("functions.jl")
include("indices.jl")
include("extract.jl")
include("interface.jl")

# TODO Add show method for ClimGrid
# TODO Refactor indices to use AxisArrays
# TODO Define a constructor for concatenation of ClimGrid type
end #module
