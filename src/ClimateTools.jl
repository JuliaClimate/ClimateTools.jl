module ClimateTools

# External modules
using Reexport
@reexport using ClimateBase
using NetCDF
@reexport using NCDatasets
using Shapefile
using AxisArrays
using NaNMath
const axes = Base.axes
using ArgCheck
using DataFrames
using Interpolations
using ProgressMeter
using Polynomials
using IterTools
using Statistics
using Random
using Dates
using GeoStats
using InverseDistanceWeighting
using Extremes
using Distances
import Base.vcat
import Base.getindex
import Base.show
import Base.size
import Base.setindex!
import Base.similar
import Base.write
import Statistics.minimum
import Statistics.maximum
import Statistics.std
import Statistics.var
import Statistics.mean
import Base: +
import Base: -
import Base: *
import Base: /
import Base.findmax

# Included files
include("functions.jl")
include("indices.jl")
include("indicators.jl")
include("extract.jl")
include("interface.jl")
include("cf_conventions.jl")
include("biascorrect.jl")
include("export.jl")
include("time.jl")
include("spatial.jl")
include("analysis.jl")

export inpoly, inpolygrid, meshgrid, inpolyvec, ndgrid
export findmax, findmin
export frostdays, summerdays, icingdays, tropicalnights
export daysabove10 #, daysbelow0, degdaysabove, degdaysbelow
export customthresover, customthresunder, annualmax, annualmin
export annualmean, annualsum, prcp1
export drought_dc
export ensemble_mean, ensemble_std, ensemble_max, ensemble_min
export load, load2D
export regrid, applymask
export shapefile_coords, shapefile_coords_poly
export resample, spatialsubset
export qqmap, qqmaptf
export biascorrect_extremes
export permute_west_east
export getdim_lat, getdim_lon, getdim_tim, isdefined, extractpoly
export get_dimname
export polyfit, polyval
export @isdefined
export merge, vaporpressure, approx_surfacepressure
export wbgt, diurnaltemperature, meantemperature
export minimum, maximum, std, var, mean
export daymean, daysum
export monthmean, monthsum, temporalmean
export yearmonthdayhour
export write, findmindist


end #module
