module ClimateTools

using CFTime
using Dates
using DimensionalData
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

include("aggregate.jl")
include("autocorrelation.jl")
include("biascorrect.jl")
include("climatology.jl")
include("ensembles.jl")
include("functions.jl")
include("markov.jl")
include("power.jl")
include("plotting.jl")
include("processERA5.jl")
include("utils.jl")


export daily_fct, climato_tp, subsample, dates_builder_yearmonth, dates_builder_yearmonth_hardcode, dates_builder_yearmonthday, dates_builder_yearmonthday_hardcode, diff, cumsum, yearly_clim
export yearly_clim
export qqmap, qqmap_bulk
export ensemble_fct
export autocorrelation
export MSModel


end
