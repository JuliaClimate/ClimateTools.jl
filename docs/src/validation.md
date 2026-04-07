# Validation and Diagnostics

ClimateTools can produce corrected scenarios, regridded fields, and large catalogs of indices, but those products still need validation before they should be trusted in downstream analysis.

## Validate Every Stage

In practice, it helps to validate at three stages:

1. after dataset preparation and subsetting
2. after regridding
3. after bias correction

## Regridding Checks

After regridding, confirm that:

- the output grid is the one you intended
- coastlines and masks look reasonable
- missing-value behavior is acceptable
- regional means or representative slices remain physically plausible

Useful quick checks:

- map a single time slice before and after regridding
- compare a regional average on the source and target grid

## Bias-Correction Checks

After bias correction, compare the corrected product with the observational reference over the calibration period.

Useful diagnostics include:

- mean bias
- standard deviation
- annual maxima and minima
- wet-day frequency or dry-spell duration
- representative grid-cell time series

Examples:

```julia
obs_txx = tx_max(obs)
qq_txx = tx_max(qq)

obs_wet = wetdays(obs_pr; thresh=1.0)
qq_wet = wetdays(qq_pr; thresh=1.0)
```

## TVC-Specific Checks

When using `tvc`, inspect not only marginal behavior but also variability and persistence.

Recommended checks:

- variance at the calibration scale
- annual maxima or heat-wave-related indices
- lag behavior or persistence-sensitive metrics at representative points

## Missing-Value Diagnostics

Many climate-processing failures show up as unexpected `NaN` propagation.

Check:

- whether the calibration domain contains missing values in key regions
- whether regridding has increased masked areas
- whether a bias-correction method returned all-`NaN` values for some cells because the calibration sample was insufficient

## Validation Is Workflow-Dependent

The right diagnostic depends on the application.

- For mean climate, annual means and standard deviations may be enough.
- For extremes, focus on maxima, wet-day intensity, spell duration, and upper-tail behavior.
- For persistence-sensitive studies, include TVC-style variability checks and duration indices.

## Related Pages

- [Bias Correction](biascorrection.md)
- [Building Climate Scenarios](scenarios.md)
- [Examples](examples.md)