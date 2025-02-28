module ClimateTools

using CircularArrays
using CFTime
using Dates
using DimensionalData
using Extremes
using LombScargle
using LongMemory
using Statistics
using YAXArrays
using NetCDF
using Zarr
using Polynomials
using Interpolations
using DataFrames
using MarSwitching
using Shapefile

include("aggregate.jl")
include("autocorrelation.jl")
include("biascorrect.jl")
include("climatology.jl")
include("ensembles.jl")
include("functions.jl")
include("gev.jl")
include("markov.jl")
include("power.jl")
include("plotting.jl")
include("processERA5.jl")
include("regrid.jl")
include("spatial.jl")
include("statistics.jl")
include("utils.jl")


export daily_fct, climato_tp, subsample, dates_builder_yearmonth, dates_builder_yearmonth_hardcode, dates_builder_yearmonthday, dates_builder_yearmonthday_hardcode, diff, cumsum, yearly_clim
export ERA5Land_dailysum
export m2mm
export yearly_clim
export quantiles
export regrid_cube
export qqmap, qqmap_bulk
export ensemble_fct
export autocorrelation
export hurst
export MSModel
export rlevels_cube


end
