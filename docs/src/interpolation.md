# Interpolation

A typical step in climate analysis is to interpolate a given grid onto another grid. `ClimateTools` provides such a tool by wrapping Scipy griddata function. It is intended for visualization or as a 1st step before bias-correcting the `ClimGrid` dataset.

[`griddata`](@ref) function will interpolate the data contained in `ClimGrid A` into the coordinates of `ClimGrid B` and returns a new `ClimGrid C` which contains the interpolated data of `A` into the grid of `B`.

```julia
C = griddata(A::ClimGrid, B::ClimGrid)
```

It is also possible to interpolate a `ClimGrid` onto specified longitude and latitude vectors and arrays.

```julia
C = griddata(A::ClimGrid, lon::AbstractArray{N, T} where N where T, lat::AbstractArray{N, T} where N where T; dimx=[], dimy=[], method::String="linear", min=[], max=[])
```

In the case a longitude and latitude 2D array is provided, the user needs to provide the dimension vectors for `x` and `y`.

## Experimental

DEPRECATED. Users will have to define their own regridding functions using GeoStats.

ClimateTools also provide a way to uses `GeoStats` geostatistics methods. See function `regrid` for some details.

```julia
using GeoStats
using ClimateTools

target = :pr # e.g. precipitation
n = 30 # max number of neighboring points
solver = Kriging(target => (maxneighbors=n,))

C = regrid(A::ClimGrid, B::ClimGrid, solver=solver)
```
