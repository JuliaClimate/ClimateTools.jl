# Home

## Overview

ClimateTools.jl is a collection of commonly-used tools in Climate science. Basics of climate field analysis will be covered, with some forays into exploratory techniques. The package is aimed to ease the typical steps of analysis climate models outputs from netCDF files that follows [Climate Forecast conventions](http://cfconventions.org/) and the creation of [climate scenarios](https://www.ouranos.ca/publication-scientifique/Guidebook-2016.pdf).

The package is registered on METADATA.jl and can be added with
```julia
] add ClimateTools
using ClimateTools
```

## Philosophy

The idea behind ClimateTools is that most, if not all, climate fields can be represented by a 2D (e.g. topography), 3D (e.g. air temperature) or 4D (e.g. winds at multiple levels) grids that are georeferenced. Those grids are named [`ClimGrid`](@ref) in ClimateTools. Every functions acts on such structure and returns a similar structure. The `ClimGrid` structure contains all elements needed to be manipulated: latitude, longitude, calendars, variable attributes, etc. that was either available in the original netCDF file or that was inferred by the metadata. Note that a `ClimGrid` is defined for a single variable.

The metadata follows the various transformations and is modified when necessary. For example, calculating the annual number of days with precipitation higher than 1mm will modify the variable name from `pr` (for precipitation) to `prcp1`, the name of the indicator. It will not, however, modify the base variable type (it will remain `pr`).

## Notes

Where possible, functions are coded to use **multiple threads**. To gain maximum performance, use (bash shell) `export JULIA_NUM_THREADS=n`, where _n_ is the number of threads. To get an idea of the number of threads you can use type (in Julia) `Sys.THREADS`. This is especially useful for climate indices, bias correction and regridding.

## Features

* Climate scenarios creation
* Extraction and visualization of CF-compliant netCDF datasets
* Custom user-provided polygons and start and end date for localized studies
* Climate indices from The joint CCl/CLIVAR/JCOMM Expert Team (ET) on Climate Change Detection and Indices (ETCCDI) as well as custom climate indices
* Regridding of a datasets onto another grid
* Post-processing of climate timeseries using Quantile-Quantile mapping method (cf. Themeßl et al. 2012, Piani et al. 2010)
* Exportation of results to a CF-compliant netCDF file
* Support for typical climate models calendars: 360_day, 365_day, Standard, Prolectip Gregorian through [NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl).
* Support for physical units through the [Unitful.jl](https://github.com/ajkeller34/Unitful.jl) package.


## Contributors

If you'd like to have other climate indices coded, please, submit them through a Pull Request! I'd be more than happy to include them. Alternatively, provide the equation in Issues.

## TO-DO

* Dashboard tool. This will return the main characteristics of a ClimGrid: maps of minimum, maximum and mean climatological values, seasonal cycle, timeseries of annual maximum, minimum and mean values, etc...
* Extreme value theory analysis

## Documentation

This documentation was built using [Documenter.jl](https://github.com/JuliaDocs).

```@example
using Dates # hide
println("Documentation built $(Dates.now()) with Julia $(VERSION)") # hide
```
