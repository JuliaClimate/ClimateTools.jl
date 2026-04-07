# Visualization

ClimateTools returns `YAXArray` and `Cube` values, so visualization is usually a matter of converting a slice to an ordinary Julia array and plotting it with your preferred graphics library.

## Quick Inspection with Plots.jl

```julia
using Plots

arr = Array(annualmax(cube))
heatmap(arr[:, :, 1], title="Annual maximum")
```

This is often enough for sanity checks after regridding, bias correction, or index calculation.

## Time-Series Inspection

ClimateTools workflows often need pointwise or regional time-series inspection in addition to maps.

```julia
series = vec(Array(cube[longitude=1, latitude=1, :]))
plot(series, title="Grid-cell time series")
```

## What to Plot During Scenario Workflows

Useful visual checks include:

- the reference observation field
- the raw model field on the same date or period
- the regridded model field
- the corrected field
- the difference between corrected and raw simulations

For bias-correction work, plotting a few representative time series is often as informative as map slices.

## Export Before Plotting

For more advanced cartography or GIS workflows, you may want to export the corrected scenario or derived index first, then use external plotting tools.

ClimateTools preserves axes and metadata as far as possible through its processing steps, so outputs remain usable as gridded products.

## Related Pages

- [Examples](examples.md)
- [Validation and Diagnostics](validation.md)
- [Building Climate Scenarios](scenarios.md)
