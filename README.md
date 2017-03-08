# ClimateTools for Julia

| **Package Status** | **Package Evaluator** | **Build Status**  |
|:------------------:|:---------------------:|:-----------------:|
| [![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip) [![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md) | [![Coverage Status](https://coveralls.io/repos/github/Balinus/ClimateTools.jl/badge.svg?branch=master)](https://coveralls.io/github/Balinus/ClimateTools.jl?branch=master) [![codecov.io](http://codecov.io/github/Balinus/ClimateTools.jl/coverage.svg?branch=master)](http://codecov.io/github/Balinus/ClimateTools.jl?branch=master) | [![Build Status](https://travis-ci.org/Balinus/ClimateTools.jl.svg?branch=master)](https://travis-ci.org/Balinus/ClimateTools.jl) [![Build status](https://ci.appveyor.com/api/projects/status/90lpp8k6430766vx?svg=true)](https://ci.appveyor.com/project/Balinus/climatetools-jl) |


This package is a collection of commonly-used tools in Climate Science. This is mainly a work-in-progress package, developed for myself and is available here, for _common-good_ purpose as well as for archive purpose. Nothing fancy here, basics of climate field analysis will be covered, with some (planned) forays into some _"state-of-the-art"_ techniques.

This package is now registered on METADATA.jl and can be added with `Pkg.add("ClimateTools")` and used with `using ClimateTools`.

The climate indices are coded to use multiple threads. To gain maximum performance, use (bash shell) `export JULIA_NUM_THREADS=n`, where _n_ is the number of threads. To get an idea of the number of threads you can use type (in Julia) `Sys.CPU_CORES`.

Since the package is evolving "rapidly", you might prefer to checkout the git repo directly, although the master might not be working.

```julia
Pkg.add("ClimateTools")
Pkg.checkout("ClimateTools")
```

## Documentation

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://balinus.github.io/ClimateTools.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://balinus.github.io/ClimateTools.jl/latest)

## Objectives

* Visualization of NetCDF files (e.g. temporal mean of a given NetCDF file), for rapid evaluation of NetCDF files
* Migration of NetCDF files to a Julia ClimGrid type
* Climate indices from The joint CCl/CLIVAR/JCOMM Expert Team (ET) on Climate Change Detection and Indices (ETCCDI)
* Custom climate indices
* Post-processing of climate timeseries using Quantile-Quantile mapping methods (cf. Piani et al. 2010)

## Examples

`using ClimateTools`

### Reading a NetCDF file
```julia
C = nc2julia(filename::String, var::String, polygon::Vector)
```

`nc2julia` return a `ClimGrid` type.

````julia
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

You can map this `ClimGrid` variable by using `mapit`:
```julia
mapit(C::ClimGrid)
```

Which should return

<p align="center">
  <img src="https://cloud.githubusercontent.com/assets/3630311/23707742/6dec18ec-03e1-11e7-90cf-0ebfd5633083.png?raw=true" alt="Precipitation example"/>
</p>

Note that if the `ClimGrid` structure has 3 dimensions (time x latitude x longitude) the `mapit` function makes a time-average.

## TO-DO

* Add a "simple" quantile-quantile mapping technique
* Add a more complex quantile-quantile mapping technique, combining POT and quantile-quantile standard technique
* Add GRIB file support (probably through [GMT.jl](https://github.com/joa-quim/GMT.jl))
