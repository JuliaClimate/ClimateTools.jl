# Indices in ClimateTools package (WIP)

## Indices

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

## Annual indices

```@docs
annualmax
```


## Functions - Tools

```@docs
ClimateTools.inpoly(p, poly::Matrix)
ClimateTools.meshgrid(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
ClimateTools.windnr(p, poly::Matrix)
ClimateTools.boxcar3(A::AbstractArray)
```

## Seasonal indices

Coming soon.
