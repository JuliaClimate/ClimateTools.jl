# Interpolation

A typical step in climate analysis is to interpolate a given grid onto another grid. `ClimateTools` provides such a tool by wrapping Scipy griddata function. It is intended for visualization or as a 1st step before bias-correcting the `ClimGrid` dataset.

`regrid` function will interpolate the data contained in `ClimGrid A` into the coordinates of `ClimGrid B` and returns a new `ClimGrid C` which contains the interpolated data of `A` into the grid of `B`.

```julia
C = regrid(A::ClimGrid, B::ClimGrid)
```

It is also possible to interpolate a `ClimGrid` onto specified longitude and latitude vectors.

```julia
C = regrid(A::ClimGrid, lon::AbstractArray{N, 1}, lat::AbstractArray{N, 1})
```
