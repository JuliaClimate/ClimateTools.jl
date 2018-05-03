# ClimateTools.jl documentation

```@contents
```

## Overview

This package is a collection of commonly-used tools in Climate Science. This is mainly a work-in-progress package, developed for myself and is available here, for _common-good_ purpose as well as for archive purpose. Nothing fancy here, basics of climate field analysis will be covered, with some forays into some _"state-of-the-art"_ techniques.

The climate indices are coded to use multiple threads. To gain maximum performance, use (bash shell) `export JULIA_NUM_THREADS=n`, where _n_ is the number of threads. To get an idea of the number of threads you can use type (in Julia) `Sys.CPU_CORES`. This can greatly reduce calculation time.

## Objectives

* Visualization of NetCDF files (e.g. temporal mean of a given NetCDF file), for rapid evaluation of NetCDF files
* Migration of NetCDF files to Julia matrix
* Climate indices from The joint CCl/CLIVAR/JCOMM Expert Team (ET) on Climate Change Detection and Indices (ETCCDI)
* Custom climate indices
* Post-processing of climate timeseries using Quantile-Quantile mapping methods (cf. Piani et al. 2010)

## Installation

This package is registered in `METADATA.jl` and so can be installed using `Pkg.add`

```julia
Pkg.add("ClimateTools")
```

For the latest version, checkout master branch

```julia
Pkg.checkout("ClimateTools")
```

## Functions - Climate indices

<!-- ```@docs
ClimateTools.frostdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.icingdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.annualmin(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.annualmax(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.tropicalnights(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.customthresover(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day}, thres)
ClimateTools.customthresunder(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day}, thres)
ClimateTools.prcp1(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.summerdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
``` -->

## Functions - Reading netCDF files

```@docs
ClimateTools.nc2julia(file::String, var::String, poly::Array{Float64})
```

```julia
type ClimGrid  
  data::AxisArray  
  model::String
  experiment::String
  run::String
  filename::String
  dataunits::String
  latunits::String
  lonunits::String
end
```



## Functions - Tools

```@docs
ClimateTools.inpoly(p, poly::Matrix)
ClimateTools.meshgrid(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.windnr(p, poly::Matrix)
ClimateTools.boxcar3(A::AbstractArray)
```

## Index

```@index
```
