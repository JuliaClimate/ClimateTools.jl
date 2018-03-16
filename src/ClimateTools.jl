module ClimateTools

# External modules
using NetCDF
using Shapefile
using AxisArrays
using ArgCheck
using PyCall
using PyPlot
using Interpolations
using ProgressMeter
import Base.vcat
import Base.getindex
import Base.show
import Base.size
import Base.endof
import Base.setindex!
import Base.similar
import Base: +
import Base: -
import Base: *
import Base: /


const basemap = PyNULL()
const np = PyNULL()
const mpl = PyNULL()
const cmocean = PyNULL()
# const scipy = PyNULL()
#const folium = PyNULL()

function __init__()
  copy!(mpl, pyimport_conda("matplotlib", "matplotlib"))
  #copy!(plt, pyimport_conda("matplotlib", "pyplot"))
  copy!(basemap, pyimport_conda("mpl_toolkits.basemap", "basemap"))
  copy!(np, pyimport_conda("numpy", "numpy"))
  copy!(cmocean, pyimport_conda("cmocean", "cmocean", "conda-forge"))
  # copy!(scipy, pyimport_conda("scipy.interpolate", "scipy"))
  #copy!(folium, pyimport_conda("folium", "folium", "conda-forge"))
  # joinpath(dirname(@__FILE__), '/Rpackages/')
  # R"install.packages('maps', lib = 'joinpath(dirname(@__FILE__), '/Rpackages/')', repo = 'http://cran.uk.r-project.org')"
end


# TYPES

struct ClimGrid{A <: AxisArray}
# struct ClimGrid
  data::A
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

end

function ClimGrid(data; model = "N/A", experiment = "N/A", run = "N/A", filename = "N/A", dataunits = "N/A", latunits = "N/A", lonunits = "N/A", variable = "N/A", typeofvar = "N/A", typeofcal = "N/A")

    ClimGrid(data, model, experiment, run, filename, dataunits, latunits, lonunits, variable, typeofvar, typeofcal)

end

# data = randn(3,2,2)
# d = 1:3
# axisdata = AxisArray(data, Axis{:time}(d), Axis{:lon}(1:2), Axis{:lat}(1:2))
# ClimGrid(axisdata)

# Included files
include("functions.jl")
include("indices.jl")
include("extract.jl")
include("interface.jl")
include("mapping.jl")
include("biascorrect.jl")

# Exported functions
export windnr, leftorright, inpoly, meshgrid, prcp1
export frostdays, summerdays, icingdays, tropicalnights
export customthresover, customthresunder, annualmax, annualmin
export annualmean, annualsum, nc2julia, sumleapyear, buildtimevec
export mapclimgrid, interp_climgrid, ClimGrid, inpolyvec, applymask
export shapefile_coords, timeresolution, pr_timefactor, spatialsubset
export temporalsubset, qqmap

end #module
