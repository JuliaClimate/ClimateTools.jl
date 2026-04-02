# Indices and Aggregations

All functions documented here accept `YAXArray` inputs. The xclim-style resampling helpers currently support `freq="YS"` and `freq="MS"`.

## Aggregation Helpers

- `daymean(cube; shifthour=0)`: daily mean.
- `daysum(cube; shifthour=0)`: daily sum.
- `yearly_resample(cube; fct=...)`: generic yearly reduction.
- `monthly_resample(cube; fct=...)`: generic monthly reduction.
- `annualmax(cube)`: annual maximum.
- `annualmin(cube)`: annual minimum.
- `annualmean(cube)`: annual mean.
- `annualsum(cube)`: annual sum.

## Legacy Annual Count Indices

- `prcp1(cube; threshold=1)`: annual count of days with precipitation at or above `threshold`.
- `frostdays(cube)`: annual count of days below `0`.
- `summerdays(cube; threshold=25)`: annual count of days above `threshold`.
- `icingdays(cube)`: annual count of days below `0`.
- `tropicalnights(cube; threshold=20)`: annual count of nights above `threshold`.
- `customthresover(cube, threshold)`: annual count of values above `threshold`.
- `customthresunder(cube, threshold)`: annual count of values below `threshold`.

## Derived Thermodynamic Indices

- `meantemperature(tasmin, tasmax)`: mean daily temperature from daily minima and maxima.
- `diurnaltemperature(tasmin, tasmax, alpha)`: weighted daily temperature estimate.
- `approx_surfacepressure(sealevel_pressure, orography, daily_temperature; temperature_unit=:auto)`: surface pressure approximation.
- `vaporpressure(specific_humidity, surface_pressure)`: vapor pressure from humidity and pressure.
- `vaporpressure(specific_humidity, sealevel_pressure, orography, daily_temperature; temperature_unit=:auto)`: vapor pressure with an internal surface-pressure approximation.
- `wbgt(mean_temperature, vapor_pressure; temperature_unit=:auto)`: simplified wet-bulb globe temperature.

## Xclim-Style Summary Indices

- `tg_max(tas; freq="YS")`: maximum daily mean temperature.
- `tg_mean(tas; freq="YS")`: mean daily mean temperature.
- `tg_min(tas; freq="YS")`: minimum daily mean temperature.
- `tx_max(tasmax; freq="YS")`: maximum daily maximum temperature.
- `tx_min(tasmax; freq="YS")`: minimum daily maximum temperature.
- `tn_max(tasmin; freq="YS")`: maximum daily minimum temperature.
- `tn_min(tasmin; freq="YS")`: minimum daily minimum temperature.
- `daily_temperature_range(tasmin, tasmax; freq="YS", op="mean")`: statistic of `tasmax - tasmin` for each period.
- `daily_temperature_range_variability(tasmin, tasmax; freq="YS")`: mean absolute day-to-day change in daily temperature range.
- `extreme_temperature_range(tasmin, tasmax; freq="YS")`: `tx_max - tn_min` for each period.
- `max_1day_precipitation_amount(pr; freq="YS")`: period maximum 1-day precipitation.
- `max_n_day_precipitation_amount(pr; window=1, freq="YS")`: period maximum rolling `window`-day precipitation sum.
- `daily_pr_intensity(pr; thresh=1.0, freq="YS", op=">=")`: mean precipitation on wet days only.

## Threshold and Run-Length Indices

- `frost_days(tasmin; thresh=0.0, freq="YS", op="<")`: count of frost days from daily minima.
- `tg_days_above(tas; thresh=10.0, freq="YS", op=">")`: count of days above a mean-temperature threshold.
- `tg_days_below(tas; thresh=10.0, freq="YS", op="<")`: count of days below a mean-temperature threshold.
- `tn_days_above(tasmin; thresh=20.0, freq="YS", op=">")`: count of nights above a minimum-temperature threshold.
- `tn_days_below(tasmin; thresh=-10.0, freq="YS", op="<")`: count of nights below a minimum-temperature threshold.
- `tx_days_above(tasmax; thresh=25.0, freq="YS", op=">")`: count of days above a maximum-temperature threshold.
- `tx_days_below(tasmax; thresh=25.0, freq="YS", op="<")`: count of days below a maximum-temperature threshold.
- `hot_days(tasmax; thresh=25.0, freq="YS", op=">")`: count of days above a hot-day threshold.
- `warm_day_frequency(tasmax; thresh=30.0, freq="YS", op=">")`: count of warm days.
- `warm_night_frequency(tasmin; thresh=22.0, freq="YS", op=">")`: count of warm nights.
- `ice_days(tasmax; thresh=0.0, freq="YS", op="<")`: count of freezing days.
- `dry_days(pr; thresh=0.2, freq="YS", op="<")`: count of dry days.
- `wetdays(pr; thresh=1.0, freq="YS", op=">=")`: count of wet days.
- `wetdays_prop(pr; thresh=1.0, freq="YS", op=">=")`: proportion of wet days.
- `growing_degree_days(tas; thresh=4.0, freq="YS")`: sum of daily mean temperature exceedances above the growing-degree threshold.
- `heating_degree_days(tas; thresh=17.0, freq="YS")`: sum of daily temperature deficits below the heating threshold.
- `cooling_degree_days(tas; thresh=18.0, freq="YS")`: sum of daily mean temperature exceedances above the cooling threshold.
- `maximum_consecutive_dry_days(pr; thresh=1.0, freq="YS")`: longest dry spell length.
- `maximum_consecutive_wet_days(pr; thresh=1.0, freq="YS")`: longest wet spell length.
- `maximum_consecutive_frost_days(tasmin; thresh=0.0, freq="YS")`: longest run of frost days.
- `maximum_consecutive_frost_free_days(tasmin; thresh=0.0, freq="YS")`: longest frost-free run.
- `maximum_consecutive_tx_days(tasmax; thresh=25.0, freq="YS")`: longest run of warm maximum-temperature days.
- `cold_spell_days(tas; thresh=-10.0, window=5, freq="YS", op="<")`: total number of days belonging to fixed-threshold cold spells.
- `cold_spell_frequency(tas; thresh=-10.0, window=5, freq="YS", op="<")`: number of fixed-threshold cold-spell events.
- `cold_spell_max_length(tas; thresh=-10.0, window=1, freq="YS", op="<")`: longest fixed-threshold cold spell.
- `cold_spell_total_length(tas; thresh=-10.0, window=3, freq="YS", op="<")`: total number of days in fixed-threshold cold spells.
- `hot_spell_frequency(tasmax; thresh=30.0, window=3, freq="YS", op=">")`: number of fixed-threshold hot-spell events.
- `hot_spell_max_length(tasmax; thresh=30.0, window=1, freq="YS", op=">")`: longest fixed-threshold hot spell.
- `hot_spell_total_length(tasmax; thresh=30.0, window=3, freq="YS", op=">")`: total number of days in fixed-threshold hot spells.
- `heat_wave_index(tasmax; thresh=25.0, window=5, freq="YS", op=">")`: total number of days in heat waves based on daily maxima alone.
- `tx_tn_days_above(tasmin, tasmax; thresh_tasmin=22.0, thresh_tasmax=30.0, freq="YS", op=">")`: count of days where daily minima and maxima both exceed their thresholds.
- `heat_wave_frequency(tasmin, tasmax; thresh_tasmin=22.0, thresh_tasmax=30.0, window=3, freq="YS", op=">")`: number of heat-wave events based on simultaneous min/max thresholds.
- `heat_wave_max_length(tasmin, tasmax; thresh_tasmin=22.0, thresh_tasmax=30.0, window=3, freq="YS", op=">")`: longest simultaneous min/max heat wave.
- `heat_wave_total_length(tasmin, tasmax; thresh_tasmin=22.0, thresh_tasmax=30.0, window=3, freq="YS", op=">")`: total number of days belonging to simultaneous min/max heat waves.
- `high_precip_low_temp(pr, tas; pr_thresh=0.4, tas_thresh=-0.2, freq="YS")`: count of days with wet precipitation thresholds and sub-threshold temperatures.

## Percentile-Based Temperature Indices

These functions expect a second `YAXArray` containing percentile thresholds already aligned to the same grid and time axis as the input temperature series.

- `tg90p(tas, tas_per; freq="YS", op=">")`: count of days where mean temperature exceeds the aligned percentile threshold.
- `tg10p(tas, tas_per; freq="YS", op="<")`: count of days where mean temperature is below the aligned percentile threshold.
- `tn90p(tasmin, tasmin_per; freq="YS", op=">")`: count of warm nights relative to an aligned percentile threshold.
- `tn10p(tasmin, tasmin_per; freq="YS", op="<")`: count of cold nights relative to an aligned percentile threshold.
- `tx90p(tasmax, tasmax_per; freq="YS", op=">")`: count of warm days relative to an aligned percentile threshold.
- `tx10p(tasmax, tasmax_per; freq="YS", op="<")`: count of cold days relative to an aligned percentile threshold.
- `cold_spell_duration_index(tasmin, tasmin_per; window=6, freq="YS", op="<")`: total number of days belonging to runs of at least `window` consecutive days below the aligned threshold.
- `warm_spell_duration_index(tasmax, tasmax_per; window=6, freq="YS", op=">")`: total number of days belonging to runs of at least `window` consecutive days above the aligned threshold.

## Examples

```julia
txx = tx_max(tasmax)
rx5day = max_n_day_precipitation_amount(pr; window=5)
wet_fraction = wetdays_prop(pr; thresh=5.0, freq="MS")

tx90 = tx90p(tasmax, tx90_threshold; freq="YS")
wsdi = warm_spell_duration_index(tasmax, tx90_threshold; window=6)
```
