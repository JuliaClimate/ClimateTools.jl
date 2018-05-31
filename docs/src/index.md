# Home

## Overview

ClimateTools.jl is a collection of commonly-used tools in Climate science. Basics of climate field analysis will be covered, with some forays into exploratory techniques. The package is aimed to ease the typical steps of analysis climate models outputs from netCDF files that follows [Climate Forecast conventions](http://cfconventions.org/) and the creation of [climate scenarios](https://www.ouranos.ca/publication-scientifique/Guidebook-2016.pdf).

The package is registered on METADATA.jl and can be added with `Pkg.add("ClimateTools")` and used with `using ClimateTools`.

## Notes

The climate indices are coded to use **multiple threads**. To gain maximum performance, use (bash shell) `export JULIA_NUM_THREADS=n`, where _n_ is the number of threads. To get an idea of the number of threads you can use type (in Julia) `Sys.CPU_CORES`.

## Features

* Climate scenarios creation
* Extraction and visualization of CF-compliant netCDF datasets
* Custom user-provided polygons and start and end date for localized studies
* Climate indices from The joint CCl/CLIVAR/JCOMM Expert Team (ET) on Climate Change Detection and Indices (ETCCDI) as well as custom climate indices
* Regridding of a datasets onto another grid
* Post-processing of climate timeseries using Quantile-Quantile mapping method (cf. Themeßl et al. 2012, Piani et al. 2010)

## Contributors

If you'd like to have other climate indices coded, please, submit them through a Pull Request! I'd be more than happy to include them. Alternatively, provide the equation in Issues.

## TO-DO

* Dashboard tool. This will return the main characteristics of a ClimGrid: maps of minimum, maximum and mean climatological values, seasonal cycle, timeseries of annual maximum, minimum and mean values, etc...
* Export ClimGrid to netCDF file.
* Add a more complex quantile-quantile mapping technique, combining extreme value theory and quantile-quantile standard technique
