# Interpolation and Regridding

Regridding is usually the bridge between raw datasets and a comparable climate-scenario workflow. Observations and model simulations often start on different grids, and most bias-correction methods are easiest to interpret once the fields share a common spatial support.

## When to Regrid

You should usually regrid when:

- observations and simulations are on different regular grids
- the simulation uses a rotated pole or curvilinear mesh
- you need a common analysis grid for multi-model comparison
- you want the corrected scenario on the observational grid

## Choosing a Method

ClimateTools exposes two main `Regridder` methods for general workflows:

- `"bilinear"` or `"linear"`: smoother, usually preferred for temperature-like continuous fields
- `"nearest_s2d"` or `"nearest"`: simpler nearest-neighbor transfer, useful when continuity is less important or for categorical behavior

As a rough rule:

- use bilinear for temperature and pressure-like continuous variables
- be more cautious with precipitation and thresholded variables, especially at coastlines or sharp gradients

## Reusable `Regridder` Workflow

Build the regridder once and reuse it across multiple variables or periods.

```julia
regridder = Regridder(source_grid_cube, destination_grid_cube; method="bilinear")

out1 = regrid(source_cube_1, regridder)
out2 = regridder(source_cube_2)
```

This is the preferred workflow when you will apply the same geometry transformation many times.

## Persisting a Regridder

```julia
save_regridder("regridder.bin", regridder)
regridder = load_regridder("regridder.bin")
```

Saved regridders validate the incoming source lon-lat grid before applying the weights. The saved format is intended for reuse within the same Julia and ClimateTools environment.

## Missing Values During Regridding

`skipna` and `na_thres` control how missing values affect the output.

```julia
out = regrid(source_cube, regridder; skipna=true, na_thres=0.25)
```

These keywords matter when the source grid contains masks, coastlines, or large missing regions.

## Dataset-Aware Workflow for Rotated or Curvilinear Grids

For datasets with a `grid_mapping` attribute, pass the dataset rather than a plain variable cube. This lets ClimateTools find the supporting 2D geographic coordinates and projection metadata.

```julia
using YAXArrays

ds = open_dataset("climex_file.nc")
dest = YAXArray(
    (Dim{:longitude}(dest_lon), Dim{:latitude}(dest_lat)),
    zeros(length(dest_lon), length(dest_lat)),
)

regridder = Regridder(ds, :pr, dest)
out = regrid(ds[:pr], regridder)
```

The dataset constructor automatically:

1. detects `grid_mapping` in the variable metadata
2. reads 2D `lon` and `lat` arrays when available
3. derives geographic coordinates from the rotated-pole metadata when needed

This is the correct route for many regional climate-model products.

## One-Shot Regridding

For one-off use, `regrid_cube` remains available.

```julia
out = regrid_cube(source_cube, destination_cube)
```

Or, with a dataset-aware source:

```julia
out = regrid_cube(ds, :pr, dest)
```

## Curvilinear and Explicit Array Interfaces

For already-geographic curvilinear grids:

```julia
out = regrid_curvilinear_to_regular(loncurv, latcurv, data, londest, latdest; method="linear")
```

For a 2D YAXArray field on a curvilinear grid:

```julia
out = regrid_curvilinear_to_regular(source_cube_2d, destination_grid_cube; method="linear")
```

For rotated grids with explicit arrays:

```julia
out = regrid_rotated_curvilinear_to_regular(lonrot, latrot, data, londest, latdest;
    grid_north_longitude=0.0,
    grid_north_latitude=90.0)
```

## Common Pitfalls

- Passing a rotated-grid variable cube directly to `Regridder(cube, dest)` instead of using the dataset-aware constructor
- Regridding after bias correction instead of before, which can make interpretation harder
- Ignoring missing-value behavior over coastlines or masked domains
- Comparing observational and model fields before confirming that the lon-lat conventions are compatible

## Role in Scenario Construction

In a climate-scenario workflow, regridding usually comes before bias correction:

1. choose the target grid
2. regrid the historical and future simulations
3. regrid the observation if needed
4. run the correction method on aligned fields

See [Building Climate Scenarios](scenarios.md) for an end-to-end example.
