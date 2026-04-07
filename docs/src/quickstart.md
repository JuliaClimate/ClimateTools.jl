# Quick Start

This page walks through the shortest useful ClimateTools workflow: open a dataset, inspect its dimensions, compute a summary, and prepare for scenario-building work.

## 1. Open a Dataset

ClimateTools operates on YAXArrays-native cubes.

```julia
using ClimateTools
using YAXArrays

cube = Cube(open_dataset("tasmax.nc"))
```

When a dataset contains multiple variables, select the one you need first.

```julia
ds = open_dataset("tasmax.nc")
tasmax = Cube(ds[:tasmax])
```

## 2. Inspect Axes and Variable Shape

Before computing anything, check which dimensions are present.

```julia
axes(tasmax)
size(tasmax)
```

Typical dimensions are `time`, `longitude`, and `latitude`, but some datasets use alternative names such as `Ti`, `rlon`, or `rlat`.

## 3. Select a Study Period

Use DimensionalData selectors to subset in time.

```julia
using Dates
using DimensionalData

hist = tasmax[time=DateTime(1981, 1, 1)..DateTime(2010, 12, 31)]
```

If your dataset uses a different time dimension name, inspect `axes(tasmax)` and adapt the selector accordingly.

## 4. Compute a Summary or Index

Aggregation functions work directly on cubes.

```julia
annual_hot = annualmax(hist)
```

Many xclim-style indices also operate on cubes directly.

```julia
txx = tx_max(hist)
hot = hot_days(hist; thresh=30.0, freq="YS")
```

## 5. Visualize a Slice

ClimateTools returns ordinary Julia arrays once you materialize the data.

```julia
using Plots

arr = Array(txx)
heatmap(arr[:, :, 1], title="Annual maximum tasmax")
```

## 6. Bias-Correct a Simulation

If you have observations, a historical simulation, and a future simulation on comparable grids, you can apply a basic bias-correction workflow.

```julia
obs = Cube(open_dataset("obs.nc"))
ref = Cube(open_dataset("historical.nc"))
fut = Cube(open_dataset("future.nc"))

qq = qqmap(obs, ref, fut; method="additive", detrend=true)
txx_future = tx_max(qq)
```

This is the core ClimateTools pattern: open, align, correct, then compute indicators.

## 7. Next Steps

- Read [Data and Subsetting](datasets.md) to understand coordinate handling.
- Read [Interpolation and Regridding](interpolation.md) before comparing datasets from different grids.
- Read [Bias Correction](biascorrection.md) and [Building Climate Scenarios](scenarios.md) for realistic scenario workflows.