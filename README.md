# ClimateTools

[![Build Status](https://travis-ci.org/Balinus/ClimateTools.jl.svg?branch=master)](https://travis-ci.org/Balinus/ClimateTools.jl)
[![Coverage Status](https://coveralls.io/repos/github/Balinus/ClimateTools.jl/badge.svg?branch=master)](https://coveralls.io/github/Balinus/ClimateTools.jl?branch=master)
[![codecov.io](http://codecov.io/github/Balinus/ClimateTools.jl/coverage.svg?branch=master)](http://codecov.io/github/Balinus/ClimateTools.jl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/90lpp8k6430766vx?svg=true)](https://ci.appveyor.com/project/Balinus/climatetools-jl)

This package is a collection of commonly-used tools in Climate Science. This is mainly a work-in-progress package, developed for myself and is available here, for _common-good_ purpose as well as for archive purpose. Nothing fancy here, basics of climate field analysis will be covered, with some forays into some _"state-of-the-art"_ techniques.

This package is now registered on METADATA.jl and can be added with `Pkg.add("ClimateTools")` and used with `using ClimateTools`.

The climate indices are coded to use multiple threads. To gain maximum performance, use (bash shell) `export JULIA_NUM_THREADS=n`, where _n_ is the number of threads. To get an idea of the number of threads you can use type (in Julia) `Sys.CPU_CORES`.

Since the package is evolving "rapidly", you might prefer to checkout the git repo directly.

`git checkout https://github.com/Balinus/ClimateTools.jl.git`

## Documentation

Minimal documentation can be found [here.](https://balinus.github.io/ClimateTools.jl/ "ClimateTools.jl documentation")

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://balinus.github.io/ClimateTools.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://balinus.github.io/ClimateTools.jl/latest)

## Objectives

* Visualization of NetCDF files (e.g. temporal mean of a given NetCDF file), for rapid evaluation of NetCDF files
* Migration of NetCDF files to Julia matrix
* Climate indices from The joint CCl/CLIVAR/JCOMM Expert Team (ET) on Climate Change Detection and Indices (ETCCDI)
* Custom climate indices
* Post-processing of climate timeseries using Quantile-Quantile mapping methods (cf. Pani et al. 2010)

## Examples

`using ClimateTools`

Rest is coming soon.

## TO-DO

* Add a "simple" quantile-quantile mapping technique
* Add a more complex quantile-quantile mapping technique, combining POT and quantile-quantile standard technique
* Add GRIB file support (probably through PyCall and pygrib package-> https://github.com/jswhit/pygrib)
