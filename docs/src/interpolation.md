# Interpolation and Regridding

ClimateTools provides YAXArrays-native regridding utilities.

## Reusable Regridder Workflow (xESMF-style)

Build a `Regridder` once and reuse it across multiple variables/time slices.

```julia
regridder = Regridder(source_grid_cube, destination_grid_cube; method="bilinear")

out1 = regrid(source_cube_1, regridder)
out2 = regridder(source_cube_2)  # equivalent callable form
```

Supported methods for `Regridder`:

- `"bilinear"` (alias `"linear"`)
- `"nearest_s2d"` (alias `"nearest"`)

You can also persist a regridder for lightweight reuse across sessions:

```julia
save_regridder("regridder.bin", regridder)
regridder = load_regridder("regridder.bin")
```

Saved regridders validate the incoming source lon/lat grid before applying the
weights. The binary format is intended for reuse within the same Julia and
ClimateTools environment.

`skipna` and `na_thres` can be passed when applying the regridder:

```julia
out = regrid(source_cube, regridder; skipna=true, na_thres=0.25)
```

## Rotated-Pole and Curvilinear Grids (Dataset path)

For datasets with a `grid_mapping` attribute (e.g. rotated-pole grids from
ClimEx/CRCM5), pass the **Dataset** instead of a plain cube so that the
`Regridder` can access the companion 2D `lon`/`lat` variables and the
`rotated_pole` metadata:

```julia
using YAXArrays

ds = open_dataset("climex_file.nc")
dest = YAXArray(
    (Dim{:longitude}(dest_lon), Dim{:latitude}(dest_lat)),
    zeros(length(dest_lon), length(dest_lat)),
)

# Build once from the Dataset
regridder = Regridder(ds, :pr, dest)

# Apply to any cube sharing the same rlon/rlat grid
out = regrid(ds[:pr], regridder)
```

The Dataset constructor automatically:
1. Detects `grid_mapping` in the variable metadata.
2. Reads 2D `lon`/`lat` arrays from the Dataset, or computes them via
   `rotated_to_geographic` using the `rotated_pole` variable attributes.
3. Builds IDW interpolation weights in geographic coordinates.

A one-shot convenience wrapper is also available:

```julia
out = regrid_cube(ds, :pr, dest)
```

The destination can also be provided as a `Dataset`. ClimateTools will prefer a
Dataset-level lon/lat coordinates directly; the destination variable name is
not used:

```julia
out = regrid_cube(climex, :pr, era5)
```

Passing a plain `YAXArray` with `rlon`/`rlat` dimensions or a `grid_mapping`
attribute directly to `Regridder(cube, dest)` will raise an error with a
message explaining how to use the Dataset path instead.

## Compatibility Wrapper

For one-off use, `regrid_cube` remains available and internally uses `Regridder`:

```julia
out = regrid_cube(source_cube, destination_cube)
```

## Curvilinear/Rotated Grid Regridding

For already-geographic curvilinear lon/lat arrays:

```julia
out = regrid_curvilinear_to_regular(loncurv, latcurv, data, londest, latdest; method="linear")
```

YAXArrays interface for 2D curvilinear fields:

```julia
out = regrid_curvilinear_to_regular(source_cube_2d, destination_grid_cube; method="linear")
```

For rotated grids with explicit arrays:

```julia
out = regrid_rotated_curvilinear_to_regular(lonrot, latrot, data, londest, latdest;
    grid_north_longitude=0.0,
    grid_north_latitude=90.0)
```
