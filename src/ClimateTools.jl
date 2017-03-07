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
using Conda
using PyCall
using PyPlot
# export PyPlot
# @pyimport mpl_toolkits.basemap as basemap
# @pyimport numpy as np

const basemap = PyNULL()
const np = PyNULL()


function __init__()
    copy!(basemap, pyimport_conda("mpl_toolkits.basemap", "basemap"))
    copy!(np, pyimport_conda("numpy", "numpy"))
end


# Exported functions
export windnr, leftorright, inpoly, meshgrid, boxcar3, prcp1, frostdays, summerdays, icingdays, tropicalnights, customthresover, customthresunder, annualmax, annualmin, nc2julia, sumleapyear, buildtimevec, inpolyV, shpextract, mapit, drawmap



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
  var::String

  function ClimGrid(data, model = "", experiment = "", run = "", filename = "", dataunits = "", latunits = "",lonunits = "", var = "")

    # to-do -> add some checks, permutedims if need be
    new(data, model, experiment, run, filename, dataunits, latunits, lonunits, var)

  end
end

# Included files
include("functions.jl")
include("indices.jl")
include("extract.jl")
include("interface.jl")
include("mapping.jl")

# TODO Add show method for ClimGrid
# TODO Refactor indices to use AxisArrays
# TODO Define a constructor for concatenation of ClimGrid type
end #module
