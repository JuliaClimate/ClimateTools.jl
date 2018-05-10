module ClimateTools

# External modules
using NetCDF
using Reexport
@reexport using NCDatasets
# using NCDatasets
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
import Base.minimum
import Base.maximum
import Base.std
import Base.var
import Base: +
import Base: -
import Base: *
import Base: /


const basemap = PyNULL()
const np = PyNULL()
const mpl = PyNULL()
const cmocean = PyNULL()
const scipy = PyNULL()
#const folium = PyNULL()

function __init__()
  copy!(mpl, pyimport_conda("matplotlib", "matplotlib"))
  #copy!(plt, pyimport_conda("matplotlib", "pyplot"))
  copy!(basemap, pyimport_conda("mpl_toolkits.basemap", "basemap"))
  copy!(np, pyimport_conda("numpy", "numpy"))
  copy!(cmocean, pyimport_conda("cmocean", "cmocean", "conda-forge"))
  copy!(scipy, pyimport_conda("scipy.interpolate", "scipy"))
  #copy!(folium, pyimport_conda("folium", "folium", "conda-forge"))
  # joinpath(dirname(@__FILE__), '/Rpackages/')
  # R"install.packages('maps', lib = 'joinpath(dirname(@__FILE__), '/Rpackages/')', repo = 'http://cran.uk.r-project.org')"
end

# TYPES
# TODO Less stringent grid type in AxisArrays. Conform more closely to CF conventions. Add information on grid type.
"""
    ClimGrid{A <: AxisArray}

In-memory representation of Climate Forecast netCDF files.
"""

struct ClimGrid{A <: AxisArray}
  data::A
  longrid::Array{N,T} where N where T
  latgrid::Array{N,T} where N where T
  msk::Array{N,T} where N where T
  grid_mapping::Dict # information of native grid
  dimension_dict::Dict
  model::String
  frequency::String
  experiment::String
  run::String
  project::String
  institute::String
  filename::String
  dataunits::String
  latunits::String # of the coordinate variable
  lonunits::String # of the coordinate variable
  variable::String # Type of variable (i.e. can be the same as "var", but it is changed when calculating indices)
  typeofvar::String # Variable type (e.g. tasmax, tasmin, pr)
  typeofcal::String # Calendar type
  varattribs::Dict # Variable attributes
  globalattribs::Dict # Global attributes

end

"""
    ClimGrid(data; longrid=[], latgrid=[], msk=[], grid_mapping=Dict(), dimension_dict=Dict(), model="NA", frequency="NA", experiment="NA", run="NA", project="NA", institute="NA", filename="NA", dataunits="NA", latunits="NA", lonunits="NA", variable="NA", typeofvar="NA", typeofcal="NA", varattribs=Dict(), globalattribs=Dict())

Constructor of the ClimGrid function. Data is an AxisArray. Everything else is optional, but usually needed for further processing (mapping, interpolation, etc...).
"""

function ClimGrid(data; longrid=[], latgrid=[], msk=[], grid_mapping=Dict(), dimension_dict=Dict(), model="NA", frequency="NA", experiment="NA", run="NA", project="NA", institute="NA", filename="NA", dataunits="NA", latunits="NA", lonunits="NA", variable="NA", typeofvar="NA", typeofcal="NA", varattribs=Dict(), globalattribs=Dict())

    if isempty(dimension_dict)
        dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
    end

    ClimGrid(data, longrid, latgrid, msk, grid_mapping, dimension_dict, model, frequency, experiment, run, project, institute, filename, dataunits, latunits, lonunits, variable, typeofvar, typeofcal, varattribs, globalattribs)

end

# function ClimGrid(data; climgrid::ClimGrid=C)
#
#     ClimGrid(data, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
#
# end

# Included files
include("functions.jl")
include("indices.jl")
include("extract.jl")
include("interface.jl")
include("mapping.jl")
include("biascorrect.jl")
include("heatsum_indices.jl")

# Exported functions
export windnr, leftorright, inpoly, inpolygrid, meshgrid, prcp1
export frostdays, summerdays, icingdays, tropicalnights
export daysabove10 #, daysbelow0, degdaysabove, degdaysbelow
export customthresover, customthresunder, annualmax, annualmin
export annualmean, annualsum, load, sumleapyear, buildtimevec
export mapclimgrid, regrid, ClimGrid, inpolyvec, applymask
export shapefile_coords, shapefile_coords_poly, timeresolution
export pr_timefactor, spatialsubset, temporalsubset
export qqmap, ndgrid, permute_west_east
export project_id, model_id, institute_id, experiment_id, frequency_var, runsim_id, getdim_lat, getdim_lon, isdefined
export @isdefined
export buildarray, timeindex, plot
export minimum, maximum, std, var

end #module
