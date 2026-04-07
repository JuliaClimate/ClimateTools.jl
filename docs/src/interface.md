# Arithmetic and Ensembles

ClimateTools is designed to feel like ordinary Julia array work on top of labeled gridded data.

## Arithmetic on Compatible Cubes

Use native broadcasting and arithmetic when cubes are already aligned.

```julia
anom = fut .- ref
ratio = fut ./ ref
```

Typical use cases include:

- anomaly fields
- ratios or percent change fields
- corrected-minus-raw comparison maps

As always, alignment matters. If the grids or time axes differ, regrid or subset before applying arithmetic.

## Ensemble Helpers

`ensemble_stats` summarizes along a chosen dimension.

```julia
stats = ensemble_stats(cube; dim="time")
stats2 = ensemble_fct(cube; dim="time")
```

`ensemble_fct` is an alias of `ensemble_stats`.

These functions are useful when:

- working with multi-member simulations
- summarizing a stack of realizations
- deriving ensemble-level diagnostics after bias correction

Those summaries now connect directly to the plotting layer. Two common patterns are:

```julia
fig = timeseriesplot(cube;
    selectors=(longitude=10, latitude=12),
    mode=:mean_ribbon)

stats = ensemble_stats(cube; dim="member")
fig2 = timeseriesplot(stats;
    selectors=(longitude=10, latitude=12),
    mode=:stats)
```

For geographic comparison of members, use `geomapfacet(cube; facetdim=:member, selectors=(time=1,))`.

## Where This Fits in the Workflow

Arithmetic and ensemble summaries are often the final layer after you have:

1. opened and aligned the data
2. regridded the datasets if needed
3. bias-corrected the simulation fields

At that stage, cube arithmetic becomes a compact way to compare scenarios.
