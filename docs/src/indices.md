# Climate Indices

## Indices

More than 20 climate indices are available in the package, such as the annual number of tropical nights, annual maximum and minimum, etc. You can calculate such indices simply with:

```julia
ind = annualmax(C::ClimGrid)
```

Which returns another `ClimGrid`. You can also map this `ClimGrid` with the `mapclimgrid` function and returns the climatological mean of the annual maximum (e.g. daily precipitation in the example below). A list of indices can be found in the documentation and in the `functions.jl` source code.

```julia
mapclimgrid(C)
```

![BNU-ESM](assets/BNU_AnnMax.png)


## Climate Indices

```@docs
annualmax
annualmean
annualmin
annualsum
approx_surfacepressure
customthresover
customthresunder
daysabove10
icingdays
frostdays
prcp1
summerdays
tropicalnights
vaporpressure
wbgt
```
