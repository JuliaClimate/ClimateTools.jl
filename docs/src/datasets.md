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

### Polygon Subsetting with spatialsubset

Use `spatialsubset` to crop and mask a cube with a polygon in lon/lat coordinates.
The polygon can be provided as `2xN` or `Nx2`.

```julia
using ClimateTools
using YAXArrays

cube = Cube(open_dataset("data.nc"))

# Coordinates are [lon; lat], NaN starts a polygon segment.
poly = [NaN -80.0 -72.0 -72.0 -80.0 -80.0;
        NaN  45.0  45.0  50.0  50.0  45.0]

sub = spatialsubset(cube, poly)
```

With shapefiles:

```julia
using ClimateTools

poly = ClimateTools.extractpoly("region.shp", n=1)
sub = spatialsubset(cube, poly)
```

Notes:

- `spatialsubset` returns a lazily-evaluated YAXArray with a bounding-box crop and polygon mask applied.
- The input cube is not fully loaded in memory for subsetting.
- The function supports longitude coordinates in either `-180..180` or `0..360` conventions.
- If the polygon does not overlap the grid, an error is raised.
