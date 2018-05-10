# Getting started

## Installation

### Required dependencies

In theory, launching the command `using ClimateTools` should install all required dependencies. However, sometimes it's just better to do it manually to ensure that all steps are properly done. If the installation fails when launching `ùsing ClimateTools`, here are the steps to do it manually.

```julia
ENV["PYTHON"] = "" # tells PyCall to use Julia's Conda python environement
Pkg.add("Conda")
using Conda
Conda.update()
Conda.add("matplotlib")
Conda.add("basemap")
Conda.add("scipy")
Conda.add("numpy")
Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.add("Pyplot")
```
### Installing ClimateTools.jl

```julia
Pkg.add("ClimateTools") # Tagged release
Pkg.checkout("ClimateTools") # For latest master branch
```

## Reading a NetCDF file

The entry point of `ClimateTools` is to load data with the `load` function. Optional polygon clipping feature is available. By providing such polygon, the `load` function  returns a `ClimGrid` with grid points contained in the polygon.

```julia
C = load(filename::String, var::String; poly::Array, data_units::String, start_date::Date, end_date::Date)
```

`load` return a `ClimGrid` type. The `ClimGrid` is a in-memory representation of a CF-compliant netCDF file for a single variable.

Using the optional `poly` argument, the user can provide a polygon and the returned `ClimGrid` will only contains the grid points inside the provided polygon. **The polygon provided should be in the -180, +180 longitude format. If the polygon crosses the International Date Line, the polygon should be splitted in multiple parts (i.e. multi-polygons).**

`start_date` and `end_date` can also be provided. It is useful when climate simulations file spans multiple decades/centuries and one only needs a temporal subset.

For some variable, the optional keyword argument `data_units` can be provided. For example, precipitation in climate models are usually provided as `kg/m^2/s`. By specifying `data_units = mm`, the `load` function returns accumulation at the data time resolution. Similarly, the user can provide `Celsius` as `data_units` and `load` will return `Celsius` instead of `Kelvin`.

```julia
struct ClimGrid
  data::AxisArray # Data
  longrid::AbstractArray{N,2} where N # the longitude grid
  latgrid::AbstractArray{N,2} where N # the latitude grid
  msk::Array{N, 2} where N # Data mask (NaNs and 1.0)
  grid_mapping::Dict#{String, Any} # bindings for native grid
  dimension_dict::Dict
  model::String
  frequency::String # Day, month, years
  experiment::String # Historical, RCP4.5, RCP8.5, etc.
  run::String
  project::String # CORDEX, CMIP5, etc.
  institute::String # UQAM, DMI, etc.
  filename::String # Path of the original file
  dataunits::String # Celsius, kelvin, etc.
  latunits::String # latitude coordinate unit
  lonunits::String # longitude coordinate unit
  variable::String # Type of variable (i.e. can be the same as "typeofvar", but it is changed when calculating indices)
  typeofvar::String # Variable type (e.g. tasmax, tasmin, pr)
  typeofcal::String # Calendar type
  varattribs::Dict # Variable attributes dictionary
  globalattribs::Dict # Global attributes dictionary
end
```

### Subsetting

Once the data is loaded in a `ClimGrid`, options to further subset the data are available.

There is a `spatialsubset` function which acts on `ClimGrid` type and further subset the data through a spatial subset using a provided polygon. The function returns a `ClimGrid`. **Polygons needs to be on a -180, +180 longitude coordinates, as data coordinates defaults to such grid.** For instance, global models are often on a 0-360 degrees grid.

```julia
C = spatialsubset(C::ClimGrid, poly:Array{N, 2} where N)
```

Temporal subset of the data is also possible with the `temporalsubset` function:

```julia
C = temporalsubset(C::ClimGrid, startdate::Date, enddate::Date)
```
