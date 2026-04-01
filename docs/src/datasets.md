# Datasets

ClimateTools uses YAXArrays for reading and manipulating gridded data.

## Open a Dataset

```julia
using YAXArrays

cube = Cube(open_dataset("data.nc"))
```

## Subset by Dimension

Use DimensionalData selectors directly:

```julia
using DimensionalData

sub = cube[Dim{:time}(DateTime(1980,1,1)..DateTime(2009,12,31))]
```

## Spatial Selection

Spatial operations in ClimateTools assume coordinate dimensions are available on the cube.
For irregular grids, use the regridding functions documented in interpolation.
