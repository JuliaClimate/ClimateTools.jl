# Function Overview

## Aggregation

- `daily_fct`
- `daymean`
- `daysum`
- `yearly_resample`
- `monthly_resample`
- `annualmax`, `annualmin`, `annualmean`, `annualsum`
- `climato_tp`
- `subsample`
- `yearly_clim`
- `ERA5Land_dailysum`

## Indices

The full index catalog lives on the dedicated [indices page](indices.md).

- Legacy annual indices: `prcp1`, `frostdays`, `summerdays`, `icingdays`, `tropicalnights`, `customthresover`, `customthresunder`
- Xclim-style summaries: `tg_max`, `tg_mean`, `tg_min`, `tx_max`, `tx_min`, `tn_max`, `tn_min`, `daily_temperature_range`, `daily_temperature_range_variability`, `extreme_temperature_range`, `max_1day_precipitation_amount`, `max_n_day_precipitation_amount`, `daily_pr_intensity`
- Threshold and spell indices: `frost_days`, `tg_days_above`, `tg_days_below`, `tn_days_above`, `tn_days_below`, `tx_days_above`, `tx_days_below`, `hot_days`, `warm_day_frequency`, `warm_night_frequency`, `ice_days`, `dry_days`, `wetdays`, `wetdays_prop`, `growing_degree_days`, `heating_degree_days`, `cooling_degree_days`, `maximum_consecutive_dry_days`, `maximum_consecutive_wet_days`, `maximum_consecutive_frost_days`, `maximum_consecutive_frost_free_days`, `maximum_consecutive_tx_days`, `cold_spell_days`, `cold_spell_frequency`, `cold_spell_max_length`, `cold_spell_total_length`, `hot_spell_frequency`, `hot_spell_max_length`, `hot_spell_total_length`, `heat_wave_index`, `tx_tn_days_above`, `heat_wave_frequency`, `heat_wave_max_length`, `heat_wave_total_length`, `high_precip_low_temp`
- Percentile-based indices: `tg90p`, `tg10p`, `tn90p`, `tn10p`, `tx90p`, `tx10p`, `cold_spell_duration_index`, `warm_spell_duration_index`

## Bias Correction

- `qqmap`
- `qqmap_bulk`
- `biascorrect_extremes`

## Time Helpers

- `timeresolution`
- `daymean_factor`
- `pr_timefactor`
- `dates_builder_yearmonth`
- `dates_builder_yearmonth_hardcode`
- `dates_builder_yearmonthday`
- `dates_builder_yearmonthday_hardcode`

## Regridding

- `Regridder`
- `regrid`
- `regrid_cube`
- `save_regridder`
- `load_regridder`
- `regrid_curvilinear_to_regular`
- `regrid_rotated_curvilinear_to_regular`

## Spatial

- `spatialsubset`

## Thermodynamics

- `vaporpressure`
- `approx_surfacepressure`
- `wbgt`
- `diurnaltemperature`
- `meantemperature`

## Ensemble

- `ensemble_stats`
- `ensemble_fct`

## Statistics and Extremes

- `quantiles`
- `rlevels_cube`

## Time Series and Regime Analysis

- `autocorrelation`
- `hurst`
- `MSModel`

## Utility Transforms

- `m2mm`
- `diff`
- `cumsum`
