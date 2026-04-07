# Troubleshooting

This page collects the most common workflow issues encountered when using ClimateTools.

## A Function Cannot Find the Time Dimension

Some functions expect a `time` axis, while others can also handle `Ti`.

What to check:

- inspect `axes(cube)`
- confirm that the dataset was opened as the intended variable
- make sure you did not drop the time axis during slicing

## Observations and Simulations Have Different Calendars

Bias-correction workflows often need comparable calibration periods and compatible calendars.

What to check:

- whether one dataset contains leap days and another does not
- whether the calibration period truly overlaps
- whether you are comparing the same variable and units

ClimateTools bias-correction methods handle common leap-day cases, but the input series still need a defensible temporal alignment.

## Regridder Fails on a Rotated Grid

If the source grid uses `rlon` and `rlat` dimensions or a `grid_mapping` attribute, do not treat it like a plain regular cube.

Use the dataset-aware path instead:

```julia
ds = open_dataset("rotated_model.nc")
regridder = Regridder(ds, :tasmax, target)
```

## Bias Correction Returns Mostly `NaN`

Possible reasons:

- the calibration sample has too many missing values
- the grids were not aligned before correction
- the time overlap is insufficient
- the variable was passed with incompatible units or wrong semantics

Check the raw data coverage first before changing the correction parameters.

## A Polygon Does Not Overlap the Grid

This often means the polygon and grid use different longitude conventions.

What to check:

- whether the grid uses `0 .. 360` but the polygon uses `-180 .. 180`
- whether the polygon coordinates are ordered correctly as lon-lat

## Regridding Looks Wrong Near Coasts

Check:

- whether `skipna` and `na_thres` should be adjusted
- whether bilinear interpolation is appropriate for the variable
- whether the source field has a sharp land-sea mask that nearest-neighbor would preserve better

## The Output Has the Right Shape but Unexpected Values

Sanity-check the workflow stage by stage:

1. inspect the raw source value at a representative grid point
2. inspect the regridded value at the same location
3. inspect the corrected value
4. compare annual or monthly summaries rather than only the full cube

## Still Unsure?

Reduce the problem to a tiny regional subset and a short time period first. That usually reveals whether the issue comes from:

- data loading
- coordinate alignment
- regridding
- bias correction
- interpretation of the derived index