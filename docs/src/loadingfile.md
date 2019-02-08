## Reading a NetCDF file

The entry point of `ClimateTools` is to load data with the `load` function. The return structure of the `load` function is a in-memory representation of the variable contained in the netCDF file.

```julia
C = load(filename::String, vari::String; poly::Array, data_units::String, start_date::Tuple, end_date::Tuple, dimension::Bool=true)
```

`load` return a `ClimGrid` type. The `ClimGrid` represent a single variable. By default, the function tries to attach physical units to the data array by using the [Unitful.jl](https://github.com/ajkeller34/Unitful.jl) package. The advantage behind physical units is that one can subtract a `ClimGrid` with `Kelvin` unit with a `ClimGrid` with `Celsius` unit and get coherent results. Be warned that some operations on some units are not allowed (you cannot "add" Celsius for instance). In the event that a user wants to do some calculations without physical logic, it is possible to load the dataset without the units by specifying `dimension=false` argument.

Using the optional `poly` argument, the user can provide a polygon and the returned `ClimGrid` will only contains the grid points inside the provided polygon. **The polygon provided should be in the -180, +180 longitude format. If the polygon crosses the International Date Line, the polygon should be splitted in multiple parts (i.e. multi-polygons).**

`start_date` and `end_date` can also be provided. It is useful when climate simulations file spans multiple decades/centuries and only a temporal subset is needed. Dates should be provided as a `Tuple` of the form `(year, month, day, hour, minute, seconds)`, where only `year` is mandatory (e.g. `(2000,)` can be provided and will defaults to `(2000, 01, 01)`).

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
  timeattrib::Dict # Time attributes (e.g. days since ... )
  varattribs::Dict # Variable attributes dictionary
  globalattribs::Dict # Global attributes dictionary
end
```

## Subsetting

Once the data is loaded in a `ClimGrid` struct, options to further subset the data are available.

### Spatial

`spatialsubset` function acts on `ClimGrid` type and subset the data through a spatial subset using a provided polygon. The function returns a `ClimGrid`. **Polygons needs to be on a -180, +180 longitude coordinates, as data coordinates defaults to such grid.** For instance, global models are often on a 0-360 degrees grid but the load function shift the data onto a -180,+180 coordinates.

```julia
C = spatialsubset(C::ClimGrid, poly:Array{N, 2} where N)
```

### Temporal

Temporal subset of the data is also possible with the `temporalsubset` function:

```julia
C = temporalsubset(C::ClimGrid, startdate::Tuple, enddate::Tuple)
```

### Discontinuous temporal (e.g. resampling)

It is also possible to only keep a given non-continuous period for a given timeframe. For example, we might be interested in keeping only northern summer months (June-July-August) from a continuous ClimGrid covering 1961-2100. `resample` returns such a subsetted ClimGrid.

```julia
Csub = resample(C, "JJA") # hardcoded ClimateTools's season
Csub = resample(C, 6, 8) # custom subset example for June-July-August
Csub = resample(C, 1, 2) # custom subset example for January-February
```
