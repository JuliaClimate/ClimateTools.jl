module ClimateTools

# External modules
using NetCDF
using Shapefile
using AxisArrays
import Base.vcat
import Base.getindex
import Base.show
import Base.size
import Base.endof
import Base.setindex!
import Base.similar
using ArgCheck
# using Conda
using PyCall
using PyPlot
# export PyPlot
# @pyimport mpl_toolkits.basemap as basemap
# @pyimport numpy as np

const basemap = PyNULL()
const np = PyNULL()
const mpl = PyNULL()

function __init__()
  copy!(mpl, pyimport_conda("matplotlib", "matplotlib"))
  copy!(basemap, pyimport_conda("mpl_toolkits.basemap", "basemap"))
  copy!(np, pyimport_conda("numpy", "numpy"))
end


# Exported functions
export windnr, leftorright, inpoly, meshgrid, prcp1, frostdays, summerdays, icingdays, tropicalnights, customthresover, customthresunder, annualmax, annualmin, annualmean, annualsum, nc2julia, sumleapyear, buildtimevec, mapclimgrid

# TYPES
immutable ClimGrid
  data::AxisArray
  model::String
  experiment::String
  run::String
  filename::String
  dataunits::String
  latunits::String
  lonunits::String
  var::String # Type of variable (i.e. can be the same as "var", but it is changed when calculating indices)
  typeofvar::String # Variable type (e.g. tasmax, tasmin, pr)
  typeofcal::String # Calendar type

  function ClimGrid(data; model = "", experiment = "", run = "", filename = "", dataunits = "", latunits = "", lonunits = "", var = "", typeofvar = "", typeofcal = "")

    # to-do -> add some checks, permutedims if need be
    new(data, model, experiment, run, filename, dataunits, latunits, lonunits, var, typeofvar, typeofcal)

  end
end

# Included files
include("functions.jl")
include("indices.jl")
include("extract.jl")
include("interface.jl")
include("mapping.jl")

end #module
