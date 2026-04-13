# Indices and Aggregations

ClimateTools provides both classic annual summaries and a broader xclim-style catalog of indices. This page groups them by analytical purpose instead of presenting a single flat list.

Unless otherwise noted, the functions documented here operate directly on `YAXArray` inputs. Many xclim-style period reductions support `freq="YS"` for yearly output and `freq="MS"` for monthly output.

## 1. Core Aggregations

Use these functions when you need a first-level summary before computing more specialized indices.

- `daymean(cube; shifthour=0)`: daily mean
- `daysum(cube; shifthour=0)`: daily sum
- `yearly_resample(cube; fct=...)`: generic yearly reduction
- `monthly_resample(cube; fct=...)`: generic monthly reduction
- `annualmax(cube)`: annual maximum
- `annualmin(cube)`: annual minimum
- `annualmean(cube)`: annual mean
- `annualsum(cube)`: annual sum

Example:

```julia
pr_day = daysum(pr)
pr_annual = annualsum(pr_day)
```

## 2. Legacy Annual Count Indices

These functions provide compact annual counts often used in older ClimateTools workflows.

- `prcp1(cube; threshold=1)`
- `frostdays(cube)`
- `summerdays(cube; threshold=25)`
- `icingdays(cube)`
- `tropicalnights(cube; threshold=20)`
- `customthresover(cube, threshold)`
- `customthresunder(cube, threshold)`

These are useful when you want a simple annual count without moving to the broader xclim-style interface.

## 3. Temperature and Precipitation Summary Indices

These functions summarize the distribution of daily values over each output period.

- `tg_max`, `tg_mean`, `tg_min`
- `tx_max`, `tx_min`
- `tn_max`, `tn_min`
- `daily_temperature_range`
- `daily_temperature_range_variability`
- `extreme_temperature_range`
- `max_1day_precipitation_amount`
- `max_n_day_precipitation_amount`
- `daily_pr_intensity`

Examples:

```julia
txx = tx_max(tasmax)
rx5day = max_n_day_precipitation_amount(pr; window=5)
wet_intensity = daily_pr_intensity(pr; thresh=1.0, freq="MS")
```

## 4. Threshold Counts

Threshold indices answer questions such as “how many hot days occurred this year?” or “how many wet days occurred each month?”

Temperature-related threshold counts:

- `frost_days`
- `tg_days_above`, `tg_days_below`
- `tn_days_above`, `tn_days_below`
- `tx_days_above`, `tx_days_below`
- `hot_days`
- `warm_day_frequency`
- `warm_night_frequency`
- `ice_days`

Precipitation-related threshold counts:

- `dry_days`
- `wetdays`
- `wetdays_prop`

Degree-day style accumulations:

- `growing_degree_days`
- `heating_degree_days`
- `cooling_degree_days`

## 5. Spell and Run-Length Indices

These functions are useful when persistence matters as much as the count itself.

- `maximum_consecutive_dry_days`
- `maximum_consecutive_wet_days`
- `maximum_consecutive_frost_days`
- `maximum_consecutive_frost_free_days`
- `maximum_consecutive_tx_days`
- `cold_spell_days`, `cold_spell_frequency`, `cold_spell_max_length`, `cold_spell_total_length`
- `hot_spell_frequency`, `hot_spell_max_length`, `hot_spell_total_length`
- `heat_wave_index`
- `tx_tn_days_above`
- `heat_wave_frequency`, `heat_wave_max_length`, `heat_wave_total_length`
- `high_precip_low_temp`

Use these when you need to quantify duration, not just frequency.

## 6. Percentile-Based Indices

These functions compare a series against a precomputed percentile threshold field.

The threshold cube must already be aligned to the same grid and time axis as the input series.

- `tg90p`, `tg10p`
- `tn90p`, `tn10p`
- `tx90p`, `tx10p`
- `cold_spell_duration_index`
- `warm_spell_duration_index`

Example:

```julia
tx90 = tx90p(tasmax, tx90_threshold; freq="YS")
wsdi = warm_spell_duration_index(tasmax, tx90_threshold; window=6)
```

## 7. Standardized Drought Indices

These functions provide a first gamma-based standardized-index workflow for monthly drought diagnostics.

- `spi`: gamma-based standardized precipitation index with zero inflation
- `spei`: gamma-first standardized precipitation evapotranspiration index for positive precomputed water-budget input

Current scope notes:

- only `freq="MS"` is supported
- input is aggregated to monthly sums before fitting
- `window` controls the rolling accumulation period in months
- `spei` currently expects a positive water-budget series; PET computation and signed-input log-logistic fitting are not part of this first pass

Example:

```julia
spi3 = spi(pr; window=3, cal_start=DateTime(1981, 1, 1), cal_end=DateTime(2010, 12, 31))
spei6 = spei(water_budget; window=6, cal_start=DateTime(1981, 1, 1), cal_end=DateTime(2010, 12, 31))
```

## 8. Thermodynamic Helpers

These are not indices in the narrow ETCCDI sense, but they are often part of scenario post-processing.

- `meantemperature(tasmin, tasmax)`
- `diurnaltemperature(tasmin, tasmax, alpha)`
- `approx_surfacepressure(sealevel_pressure, orography, daily_temperature; temperature_unit=:auto)`
- `vaporpressure(specific_humidity, surface_pressure)`
- `vaporpressure(specific_humidity, sealevel_pressure, orography, daily_temperature; temperature_unit=:auto)`
- `wbgt(mean_temperature, vapor_pressure; temperature_unit=:auto)`

## Choosing the Right Output Frequency

Many xclim-style indices allow:

- `freq="YS"` for yearly values
- `freq="MS"` for monthly values

Use yearly output for climatological summaries and monthly output when seasonal evolution matters.

## Typical Usage Pattern

In practice, indices are often computed after bias correction or on a derived scenario field.

```julia
qq = qqmap(obs, ref, fut; method="additive")
txx = tx_max(qq)
hot = hot_days(qq; thresh=30.0, freq="YS")
```

See [Examples](examples.md) and [Building Climate Scenarios](scenarios.md) for complete workflows.
