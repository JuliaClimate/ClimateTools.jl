# Interpolation and Regridding

ClimateTools provides YAXArrays-native regridding utilities.

## Regular Grid Regridding

```julia
out = regrid_cube(source_cube, destination_cube)
```

## Curvilinear/Rotated Grid Regridding

```julia
out = regrid_curvilinear_to_regular(source, dest;
    grid_north_longitude=0.0,
    grid_north_latitude=90.0)
```

For rotated grids with explicit arrays:

```julia
out = regrid_rotated_curvilinear_to_regular(lonrot, latrot, data, londest, latdest;
    grid_north_longitude=0.0,
    grid_north_latitude=90.0)
```
