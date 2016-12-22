# ClimateTools.jl documentation

```@contents
```

## Installation

This package is registered in `METADATA.jl` and so can be installed using `Pkg.add`

```julia
Pkg.add("ClimateTools")
```

For the latest version, checkout master branch

```julia
Pkg.checkout("ClimateTools")
```

## Functions

```@docs
ClimateTools.frostdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.icingdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.annualmin(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.annualmax(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.tropicalnights(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.customthresover(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day}, thres)
ClimateTools.customthresunder(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day}, thres)
ClimateTools.prcp1(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.inpoly(p, poly::Matrix)
ClimateTools.summerdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.meshgrid(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.windnr(p, poly::Matrix)
ClimateTools.boxcar3(A::AbstractArray)

```

## Index

```@index
```
