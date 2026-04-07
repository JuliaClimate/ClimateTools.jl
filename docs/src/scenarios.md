# Building Climate Scenarios

This guide describes a practical ClimateTools workflow for building climate scenarios from observations and climate-model simulations.

A typical scenario workflow has four parts:

1. Choose and prepare the observational reference.
2. Prepare historical and future climate-model simulations.
3. Interpolate all datasets to a common grid.
4. Bias-correct the model time series and compute derived products.

## Scenario Inputs

At minimum, you usually need:

- an observational or reanalysis reference series (`obs`)
- a model historical series for the calibration period (`ref` or `raw_train`)
- a model future or validation series to correct (`fut` or `raw_val`)

These may come from separate files or separate variables in the same dataset.

```julia
using ClimateTools
using YAXArrays

obs = Cube(open_dataset("obs_tasmax.nc"))
hist = Cube(open_dataset("gcm_hist_tasmax.nc"))
fut = Cube(open_dataset("gcm_ssp245_tasmax.nc"))
```

## Step 1: Check Time Axes and Calendars

Before interpolation or bias correction, inspect the time coordinates.

```julia
lookup(obs, :time)
lookup(hist, :time)
lookup(fut, :time)
```

Important checks:

- Do `obs` and `hist` cover the same calibration period?
- Do the datasets use leap days consistently?
- Are calendars compatible enough for the chosen bias-correction method?

ClimateTools bias-correction functions already handle common leap-day alignment cases by dropping February 29 where needed.

## Step 2: Select a Common Spatial Domain

If you only need a regional scenario, subset first. This reduces cost and makes diagnostics easier.

```julia
using Dates
using DimensionalData

obs_sub = obs[time=DateTime(1981, 1, 1)..DateTime(2010, 12, 31)]
hist_sub = hist[time=DateTime(1981, 1, 1)..DateTime(2010, 12, 31)]
```

You can also subset spatially with bounding boxes or polygons; see [Data and Subsetting](datasets.md).

## Step 3: Interpolate to a Common Grid

Bias correction is easiest to reason about when the observational reference and the simulation use the same grid.

For regular grids, build a reusable regridder:

```julia
regridder = Regridder(hist_sub, obs_sub; method="bilinear")

hist_rg = regrid(hist_sub, regridder)
fut_rg = regrid(fut, regridder)
```

If your source model uses a rotated pole or curvilinear grid, use the dataset-aware regridding path described on [Interpolation and Regridding](interpolation.md).

## Step 4: Choose a Bias-Correction Method

ClimateTools exposes several correction methods with different objectives.

### Quantile mapping: `qqmap`

Use `qqmap` when you want a day-of-year-based distributional correction.

```julia
qq = qqmap(obs_sub, hist_rg, fut_rg;
    method="additive",
    detrend=true,
    window=15,
    rankn=50)
```

Typical use:

- temperature with `method="additive"`
- precipitation with `method="multiplicative"`

### Bulk quantile mapping: `qqmap_bulk`

Use `qqmap_bulk` when you want a single bulk correction over the whole calibration sample rather than a day-of-year windowed correction.

### Extreme-value correction: `biascorrect_extremes`

Use `biascorrect_extremes` when you want quantile mapping plus an explicit correction of the upper tail, especially for variables where extreme behavior matters for impacts analysis.

```julia
dext = biascorrect_extremes(obs_sub, hist_rg, fut_rg; detrend=false)
```

### Time Variability Correction: `tvc`

Use `tvc` when your main objective is to correct covariance across time scales while preserving the time-event sequence of the model series.

```julia
tvc_out = tvc(obs_sub, hist_rg, fut_rg;
    scales=[365, 183, 92, 46, 23, 12, 6, 3, 2],
    eig_floor=1e-8)
```

This method is especially relevant when persistence, low-frequency variability, or heatwave-type behavior matters as much as marginal distributions.

## Step 5: Validate the Corrected Scenario

Do not stop at a successful function call. A climate scenario should be checked against the observational reference.

Recommended checks:

- compare mean and standard deviation during the calibration period
- compare annual summaries such as `annualmax`, `annualmean`, or wet-day counts
- inspect a few representative grid cells or regional aggregates
- verify that regridding and bias correction did not introduce widespread missing values

```julia
obs_txx = tx_max(obs_sub)
qq_txx = tx_max(qq)
```

More guidance is on [Validation and Diagnostics](validation.md).

## Step 6: Compute Scenario Products

Once the scenario is corrected, derive the quantities needed by the application.

Examples:

```julia
txx = tx_max(qq)
heatwaves = heat_wave_frequency(tasmin_corrected, qq; thresh_tasmin=22.0, thresh_tasmax=30.0)
rx5day = max_n_day_precipitation_amount(pr_corrected; window=5)
```

## Decision Guide

Use this rough rule of thumb:

- `qqmap`: default choice for day-of-year distribution correction
- `qqmap_bulk`: simple bulk correction without seasonal grouping
- `biascorrect_extremes`: when upper-tail behavior is central
- `tvc`: when time-scale covariance and persistence matter

## Common Failure Modes

Watch for the following before trusting a scenario:

- observations and historical simulations do not overlap in time
- datasets are on different grids and were not regridded first
- variable units are inconsistent across data sources
- large NaN masks propagate into the corrected field
- rotated-grid datasets were passed as plain cubes instead of datasets

See [Troubleshooting](troubleshooting.md) for concrete examples.

## Related Pages

- [Data and Subsetting](datasets.md)
- [Interpolation and Regridding](interpolation.md)
- [Bias Correction](biascorrection.md)
- [Indices and Aggregations](indices.md)
- [Validation and Diagnostics](validation.md)
