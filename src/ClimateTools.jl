# __precompile__()

module ClimateTools

# External modules
using NetCDF
using Reexport
@reexport using NCDatasets
# using NCDatasets
using Shapefile
using AxisArrays
const axes = Base.axes
using ArgCheck
using PyCall
using PyPlot
using Interpolations
using ProgressMeter
using Polynomials
using IterTools
using Statistics
using Dates
import Base.vcat
import Base.getindex
import Base.show
import Base.size
#import Base.endof
import Base.setindex!
import Base.similar
import Statistics.minimum
import Statistics.maximum
import Statistics.std
import Statistics.var
import Statistics.mean
import Base: +
import Base: -
import Base: *
import Base: /

const basemap = PyNULL()
const mpl = PyNULL()
const cmocean = PyNULL()
const scipy = PyNULL()

function __init__()
  copy!(mpl, pyimport_conda("matplotlib", "matplotlib"))  
  copy!(basemap, pyimport_conda("mpl_toolkits.basemap", "basemap"))
  copy!(cmocean, pyimport_conda("cmocean", "cmocean", "conda-forge"))
  copy!(scipy, pyimport_conda("scipy.interpolate", "scipy"))
end

# TYPES

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
  variable::String # Type of variable 
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

    if isempty(msk)
        msk = Array{Float64}(ones((size(data, 1), size(data, 2))))

    end

    ClimGrid(data, longrid, latgrid, msk, grid_mapping, dimension_dict, model, frequency, experiment, run, project, institute, filename, dataunits, latunits, lonunits, variable, typeofvar, typeofcal, varattribs, globalattribs)

end

"""
    TransferFunction(itp::Array, method::String, detrend::Bool)

Transfer function used during quantile-quantile mapping bias correction.
"""
struct TransferFunction
    itp::Array
    method::String
    detrend::Bool
end

# Included files
include("functions.jl")
include("indices.jl")
include("indicators.jl")
include("extract.jl")
include("interface.jl")
include("mapping.jl")
include("biascorrect.jl")
include("heatsum_indices.jl")

# Exported functions
export inpoly, inpolygrid, meshgrid, prcp1
export frostdays, summerdays, icingdays, tropicalnights
export daysabove10 #, daysbelow0, degdaysabove, degdaysbelow
export customthresover, customthresunder, annualmax, annualmin, periodmean
export annualmean, annualsum, load, load2D
export mapclimgrid, regrid, ClimGrid, inpolyvec, applymask, TransferFunction
export shapefile_coords, shapefile_coords_poly
export spatialsubset, temporalsubset, periodsubset
export qqmap, qqmaptf, ndgrid, permute_west_east
export getdim_lat, getdim_lon, isdefined, extractpoly
export polyfit, polyval
export @isdefined
export plot, merge, vaporpressure, approx_surfacepressure, wbgt, diurnaltemperature, meantemperature
export minimum, maximum, std, var, mean
export get_timevec, ensemble_mean, daymean

end #module
