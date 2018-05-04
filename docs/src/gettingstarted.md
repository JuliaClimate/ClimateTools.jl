# Getting started

## Reading a NetCDF file

The entry point of `ClimateTools` is to load data with the `nc2julia` function. Optional polygon clipping feature is available. By providing such polygon, the `nc2julia` function  returns a `ClimGrid` with grid points contained in the polygon.

```julia
C = nc2julia(filename::String, var::String; poly::Array, data_units::String, start_date::Date, end_date::Date)
```

`nc2julia` return a `ClimGrid` type. Using the optional `poly` argument, the user can provide a polygon and the returned `ClimGrid` will only contains the grid points inside the provided polygon. For some variable, the optional keyword argument `data_units` can be provided. For example, precipitation in climate models are usually provided as `kg/m^2/s`. By specifying `data_units = mm`, the `nc2julia` function returns accumulation at the data time resolution. Similarly, the user can provide `Celsius` as `data_units` and `nc2julia` will return `Celsius` instead of `Kelvin`.

The `ClimGrid` is a in-memory representation of a CF-compliant netCDF file for a single variable.

```julia
struct ClimGrid{A <: AxisArray}
# struct ClimGrid
  data::A
  longrid::AbstractArray{N,2} where N # the longitude grid
  latgrid::AbstractArray{N,2} where N # the latitude grid
  msk::Array{N, 2} where N
  grid_mapping::Dict#{String, Any} # bindings for native grid
  dimension_dict::Dict
  model::String
  frequency::String
  experiment::String
  run::String
  project::String # CORDEX, CMIP5, etc.
  institute::String
  filename::String
  dataunits::String
  latunits::String # of the coordinate variable
  lonunits::String # of the coordinate variable
  variable::String # Type of variable (i.e. can be the same as "var", but it is changed when calculating indices)
  typeofvar::String # Variable type (e.g. tasmax, tasmin, pr)
  typeofcal::String # Calendar type
  varattribs::Dict # Variable attributes
  globalattribs::Dict # Global attributes

end
```

There is a `spatialsubset` function which acts on `ClimGrid` type and further subset the data through a spatial subset using a provided polygon. The function returns a `ClimGrid`. Polygons needs to be on a -180, +180 longitude coordinates.

```julia
C = spatialsubset(C::ClimGrid, poly:Array{N, 2} where N)
```

Temporal subset of the data is also possible with the `temporalsubset` function:

```julia
C = temporalsubset(C::ClimGrid, startdate::Date, enddate::Date)
```


<!--
## Index

```@index
``` -->
