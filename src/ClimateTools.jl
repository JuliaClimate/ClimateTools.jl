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
# import Base: +
using ArgCheck
using PyCall
using PyPlot
# using Interpolations

const basemap = PyNULL()
const np = PyNULL()
const mpl = PyNULL()
# const scipy = PyNULL()
#const folium = PyNULL()

function __init__()
  copy!(mpl, pyimport_conda("matplotlib", "matplotlib"))
  #copy!(plt, pyimport_conda("matplotlib", "pyplot"))
  copy!(basemap, pyimport_conda("mpl_toolkits.basemap", "basemap"))
  copy!(np, pyimport_conda("numpy", "numpy"))
  # copy!(scipy, pyimp/ort_conda("scipy.interpolate", "scipy"))
  #copy!(folium, pyimport_conda("folium", "folium", "conda-forge"))
  # joinpath(dirname(@__FILE__), '/Rpackages/')
  # R"install.packages('maps', lib = 'joinpath(dirname(@__FILE__), '/Rpackages/')', repo = 'http://cran.uk.r-project.org')"
end


# Exported functions
export windnr, leftorright, inpoly, meshgrid, prcp1, frostdays, summerdays, icingdays, tropicalnights, customthresover, customthresunder, annualmax, annualmin, annualmean, annualsum, nc2julia, sumleapyear, buildtimevec, mapclimgrid, interp_climgrid

# TYPES
struct ClimGrid
  data::AxisArray
  model::String
  experiment::String
  run::String
  filename::String
  dataunits::String
  latunits::String
  lonunits::String
  variable::String # Type of variable (i.e. can be the same as "var", but it is changed when calculating indices)
  typeofvar::String # Variable type (e.g. tasmax, tasmin, pr)
  typeofcal::String # Calendar type

  function ClimGrid(data; model = "N/A", experiment = "N/A", run = "N/A", filename = "N/A", dataunits = "N/A", latunits = "N/A", lonunits = "N/A", variable = "N/A", typeofvar = "N/A", typeofcal = "N/A")

    # to-do -> add some checks, permutedims if need be
    new(data, model, experiment, run, filename, dataunits, latunits, lonunits, variable, typeofvar, typeofcal)

  end
end

# Included files
include("functions.jl")
include("indices.jl")
include("extract.jl")
include("interface.jl")
include("mapping.jl")

end #module
