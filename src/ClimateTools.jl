module ClimateTools

using CircularArrays
using CFTime
using Dates
using DimensionalData
using Extremes
using LombScargle
using LongMemory
using Serialization
using Statistics
using YAXArrays
using NetCDF
using Zarr
using Polynomials
using Interpolations
using NearestNeighbors
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
include("legacy_compat.jl")
include("xclim_indices.jl")
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
export yearly_resample
export monthly_resample
export daymean, daysum
export annualmax, annualmin, annualmean, annualsum
export prcp1, frostdays, summerdays, icingdays, tropicalnights, customthresover, customthresunder
export vaporpressure, approx_surfacepressure, wbgt, diurnaltemperature, meantemperature
export timeresolution, daymean_factor, pr_timefactor
export tx_max, tx_min, tn_max, tn_min, daily_temperature_range, daily_temperature_range_variability
export extreme_temperature_range, max_1day_precipitation_amount, max_n_day_precipitation_amount, daily_pr_intensity
export tg_max, tg_mean, tg_min
export frost_days, tg_days_above, tg_days_below, tn_days_above, tn_days_below, tx_days_above, tx_days_below
export hot_days, warm_day_frequency, warm_night_frequency, ice_days
export growing_degree_days, heating_degree_days, cooling_degree_days
export dry_days, wetdays, wetdays_prop, maximum_consecutive_dry_days, maximum_consecutive_wet_days
export maximum_consecutive_frost_days, maximum_consecutive_frost_free_days, maximum_consecutive_tx_days
export cold_spell_days, cold_spell_frequency, cold_spell_max_length, cold_spell_total_length
export hot_spell_frequency, hot_spell_max_length, hot_spell_total_length, heat_wave_index
export tx_tn_days_above, heat_wave_frequency, heat_wave_max_length, heat_wave_total_length, high_precip_low_temp
export tg90p, tg10p, tn90p, tn10p, tx90p, tx10p
export cold_spell_duration_index, warm_spell_duration_index
export quantiles
export Regridder, regrid, regrid_cube, save_regridder, load_regridder
export regrid_curvilinear_to_regular, regrid_rotated_curvilinear_to_regular
export qqmap, qqmap_bulk, biascorrect_extremes
export ensemble_fct, ensemble_stats
export autocorrelation
export hurst
export MSModel
export rlevels_cube
export spatialsubset


end
