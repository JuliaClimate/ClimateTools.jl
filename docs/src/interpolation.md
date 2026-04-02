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
