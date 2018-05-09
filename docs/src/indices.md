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

Climate Indices functions also accept other type of arguments. For example, `annualmax` can be called with the following type:

```julia
ind = annualmax(data::Array{Float64, 3}, dates::StepRange{Date, Base.Dates.Day})
```

Although, without metadata contained into a `ClimGrid`, no checks are done.


## Annual indices

```@docs
annualmax
```

## Seasonal indices

Coming soon.
