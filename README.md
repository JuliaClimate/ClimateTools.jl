# ClimateTools for Julia

| **Package Status** | **Package Evaluator** | **Build Status**  |
|:------------------:|:---------------------:|:-----------------:|
| [![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) [![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md) | [![Coverage Status](https://coveralls.io/repos/github/Balinus/ClimateTools.jl/badge.svg?branch=master)](https://coveralls.io/github/Balinus/ClimateTools.jl?branch=master) [![codecov.io](http://codecov.io/github/Balinus/ClimateTools.jl/coverage.svg?branch=master)](http://codecov.io/github/Balinus/ClimateTools.jl?branch=master) | [![Build Status](https://travis-ci.org/Balinus/ClimateTools.jl.svg?branch=master)](https://travis-ci.org/Balinus/ClimateTools.jl) [![Build status](https://ci.appveyor.com/api/projects/status/90lpp8k6430766vx?svg=true)](https://ci.appveyor.com/project/Balinus/climatetools-jl) |

N.B. master branch should be compatible with Julia 0.6.

## Documentation

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://balinus.github.io/ClimateTools.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://balinus.github.io/ClimateTools.jl/latest)

This package is a collection of commonly-used tools in Climate Science. This is mainly a work-in-progress package, developed for myself and is available here, for _common-good_ purpose as well as for archive purpose. Basics of climate field analysis will be covered, with some (planned) forays into some _state-of-the-art_ techniques.

This package is registered on METADATA.jl and can be added with `Pkg.add("ClimateTools")` and used with `using ClimateTools`.

The climate indices are coded to use multiple threads. To gain maximum performance, use (bash shell) `export JULIA_NUM_THREADS=n`, where _n_ is the number of threads. To get an idea of the number of threads you can use type (in Julia) `Sys.CPU_CORES`.

_Note: using nested parralel loops will not give the expected results (not thread safe). For example, if you already have multiple threaded loops in your code and you plan to call `ClimateTools.jl` indices, this won't give the expected results. Use this package with caution in that case. In a future version, the multiple threads will be an optional parameter, giving the user the option to integrate the package with their framework._

Since the package is evolving "rapidly", you might prefer to checkout the git repo directly, although the master might not be working (I usually don't push broken version though).

```julia
Pkg.add("ClimateTools")
Pkg.checkout("ClimateTools")
```

## Contributors

If you'd like to have other climate indices coded, please, submit them through a Pull Request! I'd be more than happy to include them. Alternatively, provide the equation in Issues.

## Objectives

* Visualization of NetCDF files (e.g. temporal mean of a given NetCDF file), for rapid evaluation of NetCDF files
* Migration of NetCDF files to a Julia `ClimGrid` type
* Climate indices from The joint CCl/CLIVAR/JCOMM Expert Team (ET) on Climate Change Detection and Indices (ETCCDI)
* Custom climate indices
* Post-processing of climate timeseries using Quantile-Quantile mapping methods (cf. Piani et al. 2010)

## Examples

`using ClimateTools`

### Reading a NetCDF file
```julia
C = nc2julia(filename::String, var::String; poly::Array)
```

`nc2julia` return a `ClimGrid` type.

```julia
immutable ClimGrid
  data::AxisArray
  model::String
  experiment::String
  run::String
  filename::String
  dataunits::String
  latunits::String
  lonunits::String
  var::String
end
```

### Merging ClimGrid type
Sometimes, the timeseries are split among multiple files. To obtain the complete timeseries, you can `merge` 2 `ClimGrid`. The method is based on the merging of 2 `AxisArrays` and is overloaded for the `ClimGrid` type.

```julia
C = merge(C1::ClimGrid, C2::ClimGrid)
```

### Mapping the ClimGrid type

You can map this `ClimGrid` variable by using `mapclimgrid`:
```julia
mapclimgrid(C::ClimGrid; region = "World")
```

Which should return

<p align="center">
  <img src="https://cloud.githubusercontent.com/assets/3630311/23712122/e97bd322-03ef-11e7-93da-749c961c4070.png?raw=true" width="771" height="388" alt="Precipitation example"/>
</p>

Note that if the `ClimGrid` data structure has 3 dimensions (time x latitude x longitude) the `mapclimgrid` function makes a time-average (i.e. climatological mean). Right now, 3 options are available for region: `World`, `Canada` and the default `auto` which use the maximum and minimum of the lat-long coordinates inside the `ClimGrid` structure. Other regions will be added in the future, as well as the option to send a custom region defined by a lat-lon box.

In a future release, the user will have the option to specify his own time period (e.g. plotting the time-average of a given month and year, as opposed to the time-average of the whole `ClimGrid` structure).

### Indices

Some indices are available in the packages, such as the annual number of tropical nights, annual maximum and minimum, etc. you can calculate such indices with:

```julia
ind = annualmax(C::ClimGrid)
```

Which returns another `ClimGrid`. You can also map this `ClimGrid` with the `mapclimgrid` function and returns the climatological mean of the annual maximum (e.g. daily precipitation in the example below). A list of indices can be found in the documentation.

<p align="center">
  <img src="https://cloud.githubusercontent.com/assets/3630311/23873133/59b85c08-0807-11e7-967b-7cc7d28aada0.png?raw=true" width="771" height="388" alt="Precipitation example"/>
</p>


## TO-DO

* Add a "simple" quantile-quantile mapping technique
* Add a more complex quantile-quantile mapping technique, combining POT and quantile-quantile standard technique
* Add GRIB file support (probably through [GMT.jl](https://github.com/joa-quim/GMT.jl))
