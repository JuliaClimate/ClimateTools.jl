# Bias Correction

ClimateTools provides several bias-correction methods aimed at different goals. Some methods mainly correct distributions, while others target time-scale variability or upper-tail behavior.

Before choosing a method, make sure that:

- the observational reference and the model simulation cover a comparable calibration period
- the datasets have compatible units
- the grids are compatible, or have already been regridded to a common target

## Which Method Should You Use?

Use this quick guide:

- `qqmap`: default day-of-year quantile mapping for many temperature and precipitation workflows
- `qqmap_bulk`: bulk quantile mapping without day-of-year grouping
- `biascorrect_extremes`: upper-tail-aware correction when extremes matter
- `tvc`: correction of temporal covariance across time scales while preserving event order

## Quantile Mapping with `qqmap`

`qqmap` performs day-of-year-based quantile mapping.

```julia
qq = qqmap(obs, ref, fut;
    method="additive",
    detrend=true,
    window=15,
    rankn=50,
    qmin=0.01,
    qmax=0.99)
```

Typical choices:

- `method="additive"` for temperature-like variables
- `method="multiplicative"` for precipitation-like variables

Important parameters:

- `window`: day-of-year neighborhood used for seasonal grouping
- `rankn`: number of quantiles used in the mapping
- `detrend`: whether to remove and reapply a polynomial trend during correction

Use `qqmap` when seasonal behavior matters and you want the correction to vary through the annual cycle.

## Bulk Quantile Mapping with `qqmap_bulk`

```julia
qq_bulk = qqmap_bulk(obs, ref, fut; method="multiplicative")
```

This method uses the full sample rather than a moving day-of-year window. It is simpler, but it does not adapt the correction seasonally.

## Extreme-Value Tail Correction with `biascorrect_extremes`

Use `biascorrect_extremes` when the upper tail needs special treatment beyond standard quantile mapping.

```julia
dext = biascorrect_extremes(obs, ref, fut;
    detrend=false,
    P=0.95,
    runlength=2,
    frac=0.25,
    power=1.0)
```

The function combines multiplicative quantile mapping with an explicit generalized Pareto distribution correction of extreme values.

The function also accepts optional external extreme-value parameters with columns `lat`, `lon`, `mu`, `sigma`, and `xi`:

```julia
using DataFrames

gevparams = DataFrame(
    lat=[45.0, 46.0],
    lon=[-73.0, -72.0],
    mu=[20.0, 21.5],
    sigma=[4.0, 4.5],
    xi=[0.10, 0.08],
)

dext = biascorrect_extremes(obs, ref, fut; gevparams=gevparams)
```

This method is most appropriate when impact-relevant extremes are central to the study.

## Time Variability Correction with `tvc`

`tvc` applies the Time Variability Correction method of Shao et al. (2024), DOI 10.1029/2023MS003640.

The method is designed to correct covariance across and between multiple time scales while preserving the event sequence of the validation series.

```julia
corrected = tvc(obs, raw_train, raw_val;
    scales=[365, 183, 92, 46, 23, 12, 6, 3, 2],
    eig_floor=1e-8,
    keep_original=false)
```

For repeated applications to multiple validation series, fit once and reuse the fitted model:

```julia
model = fit_tvc(obs_series, raw_train_series; scales=[365, 183, 92, 46, 23, 12, 6, 3, 2])
corrected_series = apply_tvc(model, raw_val_series)
```

`tvc` returns a series on the validation time axis. The leading warm-up segment, of length `sum(scales) - length(scales)`, is filled with `NaN` because the multiscale rolling decomposition needs prior samples before the corrected series is well-defined.

## Practical Workflow

For scenario construction, the usual order is:

1. subset the calibration period
2. regrid observations and simulations to a common grid
3. choose a correction method
4. validate the corrected output before computing derived indicators

See [Building Climate Scenarios](scenarios.md) for the full workflow.

## Validation Checklist

After correction, compare the corrected data against the observational reference over the calibration period.

Useful checks include:

- changes in mean and standard deviation
- annual maxima and minima
- wet-day fractions or dry-spell lengths
- map patterns at representative dates
- time-series slices at representative grid points

See [Validation and Diagnostics](validation.md) for more detailed guidance.
