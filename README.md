# ClimateTools for Julia

| **Package Status** | **Package Evaluator** | **Build Status**  |  **DOI**  |
|:------------------:|:---------------------:|:-----------------:|:---------:|
| [![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)  | [![Coverage Status](https://coveralls.io/repos/github/Balinus/ClimateTools.jl/badge.svg?branch=master)](https://coveralls.io/github/Balinus/ClimateTools.jl?branch=master) [![codecov.io](http://codecov.io/github/Balinus/ClimateTools.jl/coverage.svg?branch=master)](http://codecov.io/github/Balinus/ClimateTools.jl?branch=master) | [![Build Status](https://travis-ci.org/Balinus/ClimateTools.jl.svg?branch=master)](https://travis-ci.org/Balinus/ClimateTools.jl)| [![DOI](https://zenodo.org/badge/76293821.svg)](https://zenodo.org/badge/latestdoi/76293821) |



## Documentation

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://balinus.github.io/ClimateTools.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://balinus.github.io/ClimateTools.jl/latest)

This package is a collection of commonly-used tools in Climate Science. Basics of climate field analysis will be covered, with some forays into exploratory techniques. The package is aimed to ease the typical steps of analysis climate models outputs and observed time series from weather stations.

This package is registered on METADATA.jl and can be added with `Pkg.add("ClimateTools")` and used with `using ClimateTools`.

```julia
Pkg.add("ClimateTools") # Tagged release
Pkg.checkout("ClimateTools") # For latest master branch
```

The climate indices are coded to use **multiple threads**. To gain maximum performance, use (bash shell) `export JULIA_NUM_THREADS=n`, where _n_ is the number of threads. To get an idea of the number of threads you can use type (in Julia) `Sys.CPU_CORES`.

## Features

* Extraction and visualization of CF-compliant netCDF datasets
* Custom user-provided polygons and start and end date for localized studies
* Climate indices from The joint CCl/CLIVAR/JCOMM Expert Team (ET) on Climate Change Detection and Indices (ETCCDI) as well as custom climate indices
* Regridding of a datasets onto another grid
* Post-processing of climate timeseries using Quantile-Quantile mapping method (cf. Themeßl et al. 2012, Piani et al. 2010)

## Getting started

`using ClimateTools`

### Reading a NetCDF file

The entry point of `ClimateTools` is to load data with the `nc2julia` function. Optional polygon clipping feature is available. By providing such polygon, the `nc2julia` function  returns a `ClimGrid` with grid points contained in the polygon.

```julia
C = nc2julia(filename::String, var::String; poly::Array, data_units::String, start_date::Date, end_date::Date)
```

`nc2julia` return a `ClimGrid` type. Using the optional `poly` argument, the user can provide a polygon and the returned `ClimGrid` will only contains the grid points inside the provided polygon. For some variable, the optional keyword argument `data_units` can be provided. For example, precipitation in climate models are usually provided as `kg/m^2/s`. By specifying `data_units = mm`, the `nc2julia` function returns accumulation at the data time resolution. Similarly, the user can provide `Celsius` as `data_units` and `nc2julia` will return `Celsius` instead of `Kelvin`.

The `ClimGrid` is a in-memory representation of a CF-compliant netCDF file for a single variable.

```julia
struct ClimGrid
  data::AxisArray
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

Furthermore, there is also the `spatialsubset` function which acts on `ClimGrid` type and further subset the data through a spatial subset using a user polygon. The function returns a `ClimGrid`.

```julia
C = spatialsubset(C::ClimGrid, poly:Array{N, 2} where N)
```

Temporal subset of the data is also possible with the `temporalsubset` function:

```julia
C = temporalsubset(C::ClimGrid, startdate::Date, enddate::Date)
```

### Mapping the ClimGrid type

Mapping climate information can be done by using `mapclimgrid`:
```julia
mapclimgrid(C::ClimGrid; region = "World")
```

Which should return

<p align="center">
  <img src="https://cloud.githubusercontent.com/assets/3630311/23712122/e97bd322-03ef-11e7-93da-749c961c4070.png?raw=true" width="771" height="388" alt="Precipitation example"/>
</p>

Note that if the `ClimGrid` data structure has 3 dimensions (time x longitude x latitude) the `mapclimgrid` function makes a time-average (i.e. climatological mean). Right now, options are available for region: `World`, `Canada`, `Quebec` and the default `auto` which use the maximum and minimum of the lat-long coordinates inside the `ClimGrid` structure. The user can also provide a polygon(s) and the `mapclimgrid` function will clip the grid points outside the specified polygon. Another option is to provide a mask (with dimensions identical to the spatial dimension of the `ClimGrid` data) which contains `NaN` and `1.0` and the data inside the `ClimGrid` struct will be clipped with the mask. Other regions will be added in the future, as well as the option to send a custom region defined by a lat-lon box.

In a future release, the user will have the option to specify his own time period (e.g. plotting the time-average of a given month and year, as opposed to the time-average of the whole `ClimGrid` structure).

### Indices

More than 20 climate indices are available in the package, such as the annual number of tropical nights, annual maximum and minimum, etc. You can calculate such indices simply with:

```julia
ind = annualmax(C::ClimGrid)
```

Which returns another `ClimGrid`. You can also map this `ClimGrid` with the `mapclimgrid` function and returns the climatological mean of the annual maximum (e.g. daily precipitation in the example below). A list of indices can be found in the documentation and in the `functions.jl` source code.

<p align="center">
  <img src="https://cloud.githubusercontent.com/assets/3630311/23873133/59b85c08-0807-11e7-967b-7cc7d28aada0.png?raw=true" width="771" height="388" alt="Precipitation example"/>
</p>

Climate Indices functions also accept other type of argument. For example, `annualmax` can be called with the following type:

```julia
ind = annualmax(data::Array{Float64, 3}, dates::StepRange{Date, Base.Dates.Day})
```

### Interpolation

A typical step in climate analysis is to interpolate a given grid onto another grid. `ClimateTools` provides such a tool by wrapping Scipy griddata function. It is intended for visualization or as a 1st step before bias-correcting the `ClimGrid` dataset.

The following command will interpolate the data contained in `ClimGrid A` into the coordinates of `ClimGrid B` and returns a new `ClimGrid C` which contains the interpolated data of `A` into the grid of `B`.

```julia
C = interp_climgrid(A::ClimGrid, B::ClimGrid)
```

It is also possible to interpolate a `ClimGrid` onto specified longitude and latitude vectors.

```julia
C = interp_climgrid(A::ClimGrid, lon::AbstractArray{N, 1}, lat::AbstractArray{N, 1})
```

### Bias-correction

TODO

### Merging ClimGrid type

Sometimes, the timeseries are split among multiple files (e.g. climate models outputs). To obtain the complete timeseries, you can `merge` 2 `ClimGrid`. The method is based on the merging of 2 `AxisArrays` and is overloaded for the `ClimGrid` type.

```julia
C = merge(C1::ClimGrid, C2::ClimGrid)
```

## Contributors

If you'd like to have other climate indices coded, please, submit them through a Pull Request! I'd be more than happy to include them. Alternatively, provide the equation in Issues.

## TO-DO

* Dashboard tool. This will return the main characteristics of a ClimGrid: maps of minimum, maximum and mean climatological values, seasonal cycle, timeseries of annual maximum, minimum and mean values, etc...
* Create a WeatherStation type.
* Export ClimGrid to netCDF file.
* Add a more complex quantile-quantile mapping technique, combining extreme value theory and quantile-quantile standard technique


N.B. version 0.1.2 is compatible with Julia 0.5 and version >0.2.0 is for Julia 0.6. To use a specific version of the package, you can use in Julia the following command:

```julia
Pkg.pin("ClimateTools",v"0.1.2") # if using Julia 0.5
```
