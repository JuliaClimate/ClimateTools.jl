# Data and Subsetting

ClimateTools uses YAXArrays and DimensionalData selectors for reading and manipulating gridded climate datasets.

## Opening Datasets

The usual entry point is `open_dataset`, followed by `Cube` when you want a single variable.

```julia
using YAXArrays

cube = Cube(open_dataset("data.nc"))
```

When a file contains several variables, keep the dataset object and select explicitly.

```julia
ds = open_dataset("data.nc")
tasmax = Cube(ds[:tasmax])
pr = Cube(ds[:pr])
```

This is often clearer in climate-scenario workflows where the observational and simulation files contain several candidate variables.

## Inspecting Axes

Before computing on a cube, inspect the axes and dimension names.

```julia
axes(tasmax)
size(tasmax)
```

ClimateTools functions usually expect a time dimension and one or more spatial dimensions. Common names are:

- `time` or `Ti`
- `longitude` / `latitude`
- `lon` / `lat`
- `rlon` / `rlat` for rotated grids

Many workflows fail not because of the data values themselves, but because time or spatial coordinates are missing or named unexpectedly.

## Selecting Time Ranges

Use DimensionalData selectors directly.

```julia
using Dates
using DimensionalData

hist = tasmax[time=DateTime(1980, 1, 1)..DateTime(2009, 12, 31)]
```

If the cube uses a different time dimension name, inspect `axes(cube)` and adapt the selector.

## Selecting Spatial Ranges

For regular lon-lat grids, you can subset the same way.

```julia
regional = tasmax[
    longitude=-80.0..-60.0,
    latitude=45.0..55.0,
]
```

This is usually a good first step before regridding or bias correction when you only need a study region.

## Coordinate Conventions

Climate datasets often mix longitude conventions:

- `-180 .. 180`
- `0 .. 360`

ClimateTools spatial utilities handle both conventions in many cases, but you should still verify that:

- the polygon or bounding box is expressed in the same convention as the grid
- the observational reference and model grid are geographically compatible before regridding

## Polygon Subsetting with `spatialsubset`

Use `spatialsubset` to crop and mask a cube with a polygon in lon-lat coordinates.

The polygon can be provided as `2xN` or `Nx2`, with `NaN` separating polygon segments.

```julia
using ClimateTools
using YAXArrays

cube = Cube(open_dataset("data.nc"))

poly = [NaN -80.0 -72.0 -72.0 -80.0 -80.0;
        NaN  45.0  45.0  50.0  50.0  45.0]

sub = spatialsubset(cube, poly)
```

With shapefiles:

```julia
poly = ClimateTools.extractpoly("region.shp", n=1)
sub = spatialsubset(cube, poly)
```

Notes:

- `spatialsubset` returns a lazily-evaluated `YAXArray` with a bounding-box crop and polygon mask applied.
- The input cube is not fully loaded into memory during the subsetting step.
- If the polygon does not overlap the grid, the function raises an error.

## Missing Data and Masks

ClimateTools functions usually propagate missing or `NaN` values rather than guessing replacements. Before bias correction or index calculation, inspect whether the domain has:

- permanent land-sea masks
- missing months or days
- irregular missingness in the calibration period

Large masked areas can affect regridding and bias-correction robustness.

## Rotated and Curvilinear Grids

If your simulation uses `rlon` and `rlat` dimensions or has a `grid_mapping` attribute, prefer dataset-aware regridding workflows rather than manual slicing of a plain variable cube. See [Interpolation and Regridding](interpolation.md).

## Role in Scenario Workflows

Dataset preparation is the first stage of climate-scenario construction. A good rule is:

1. open the observation and model datasets
2. inspect the dimensions and calendars
3. subset the study period and region
4. only then regrid and bias-correct

See [Building Climate Scenarios](scenarios.md) for the full workflow.
