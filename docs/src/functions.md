# API Overview

This page groups the exported ClimateTools API by workflow. Use it as a map between the narrative guides and the concrete public functions.

## Aggregation and Resampling

- `daily_fct`: generic daily reduction over the time dimension
- `daymean`, `daysum`: daily mean and daily sum
- `yearly_resample`, `monthly_resample`: generic grouped time reductions
- `annualmax`, `annualmin`, `annualmean`, `annualsum`: annual summaries
- `climato_tp`, `subsample`, `yearly_clim`, `ERA5Land_dailysum`: helper functions for specific aggregation workflows

See [Indices and Aggregations](indices.md).

## Bias Correction

- `qqmap`: day-of-year quantile mapping
- `qqmap_bulk`: bulk quantile mapping
- `biascorrect_extremes`: extreme-tail-aware bias correction
- `TVCModel`, `fit_tvc`, `apply_tvc`, `tvc`: time variability correction API

See [Bias Correction](biascorrection.md) and [Building Climate Scenarios](scenarios.md).

## Regridding and Spatial Alignment

- `Regridder`: reusable regridding object
- `regrid`: apply a saved or precomputed regridder
- `regrid_cube`: one-shot regridding wrapper
- `save_regridder`, `load_regridder`: persistence helpers
- `regrid_curvilinear_to_regular`, `regrid_rotated_curvilinear_to_regular`: explicit curvilinear and rotated-grid interfaces

See [Interpolation and Regridding](interpolation.md).

## Spatial Subsetting

- `spatialsubset`: polygon-based crop and mask
- Regular lon/lat cubes can be subset directly; rotated or curvilinear grids should be passed as a parent `Dataset`

See [Data and Subsetting](datasets.md).

## Indices

The full grouped catalog is described on [Indices and Aggregations](indices.md).

- Legacy annual indices: `prcp1`, `frostdays`, `summerdays`, `icingdays`, `tropicalnights`, `customthresover`, `customthresunder`
- Xclim-style summaries: `tg_max`, `tg_mean`, `tg_min`, `tx_max`, `tx_min`, `tn_max`, `tn_min`, `daily_temperature_range`, `daily_temperature_range_variability`, `extreme_temperature_range`, `max_1day_precipitation_amount`, `max_n_day_precipitation_amount`, `daily_pr_intensity`
- Threshold and spell indices: `frost_days`, `tg_days_above`, `tg_days_below`, `tn_days_above`, `tn_days_below`, `tx_days_above`, `tx_days_below`, `hot_days`, `warm_day_frequency`, `warm_night_frequency`, `ice_days`, `dry_days`, `wetdays`, `wetdays_prop`, `growing_degree_days`, `heating_degree_days`, `cooling_degree_days`, `maximum_consecutive_dry_days`, `maximum_consecutive_wet_days`, `maximum_consecutive_frost_days`, `maximum_consecutive_frost_free_days`, `maximum_consecutive_tx_days`, `cold_spell_days`, `cold_spell_frequency`, `cold_spell_max_length`, `cold_spell_total_length`, `hot_spell_frequency`, `hot_spell_max_length`, `hot_spell_total_length`, `heat_wave_index`, `tx_tn_days_above`, `heat_wave_frequency`, `heat_wave_max_length`, `heat_wave_total_length`, `high_precip_low_temp`
- Occurrence and season indices: `first_day_temperature_above`, `first_day_temperature_below`, `first_snowfall`, `last_snowfall`, `last_spring_frost`, `growing_season_start`, `growing_season_end`, `growing_season_length`, `frost_free_season_start`, `frost_free_season_end`, `frost_free_season_length`, `frost_season_length`
- Percentile-based indices: `tg90p`, `tg10p`, `tn90p`, `tn10p`, `tx90p`, `tx10p`, `cold_spell_duration_index`, `warm_spell_duration_index`
- Standardized drought indices: `spi`, `spei`

## Thermodynamic Helpers

- `vaporpressure`
- `approx_surfacepressure`
- `wbgt`
- `diurnaltemperature`
- `meantemperature`

These are often used after bias correction or in post-processing chains.

## Ensembles and Statistics

- `ensemble_stats`, `ensemble_fct`: ensemble summary helpers
- `ensemble_mean_std_max_min`: xclim-style ensemble mean, standard deviation, maximum, and minimum over `realization`
- `ensemble_percentiles`: percentile summaries over `realization`, with split or stacked outputs and optional weights
- `make_criteria`: flatten multi-variable ensemble diagnostics to a 2D `realization x criteria` cube for reduction workflows
- `kkz_reduce_ensemble`: deterministic KKZ ensemble-member subset selection
- `robustness_fractions`, `robustness_categories`, `robustness_coefficient`: sign-agreement and robustness diagnostics for projected change
- `quantiles`: grouped quantile summaries over a selected dimension
- `gevfit_cube`, `gpfit_cube`, `returnlevel_cube`: reusable extreme-value fitting and return-level reuse over a reduced cube dimension; `gevfit_cube` is for block-maximum inputs, while `gpfit_cube` is for threshold exceedances
- `rlevels_cube`: legacy one-step return-level estimation on raw cubes

## Plotting

- `geomap`: single-panel geographic plotting for regular and rotated-grid data
- `geomapfacet`: faceted geographic plotting over time or ensemble-like dimensions
- `timeseriesplot`: line, multi-member, ribbon, `ensemble_stats`-aware, and dataset-variable time-series plots
- `statsplot`: histogram, boxplot, scatter, grouped ensemble-distribution, and dataset-variable plots
- `robustnessmap`: categorical geographic plotting for robustness fractions or robustness-category outputs

## Time-Series and Regime Diagnostics

- `autocorrelation`
- `hurst`
- `MSModel`

## Time Helpers and Utility Transforms

- `timeresolution`
- `daymean_factor`
- `pr_timefactor`
- `dates_builder_yearmonth`
- `dates_builder_yearmonth_hardcode`
- `dates_builder_yearmonthday`
- `dates_builder_yearmonthday_hardcode`
- `m2mm`
- `diff`
- `cumsum`

## Suggested Reading Order

If you arrived here before reading the guides, start with:

1. [Quick Start](quickstart.md)
2. [Data and Subsetting](datasets.md)
3. [Interpolation and Regridding](interpolation.md)
4. [Bias Correction](biascorrection.md)
5. [Building Climate Scenarios](scenarios.md)
