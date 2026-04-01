# Interface

ClimateTools functions are defined on YAXArray/Cube values.

## Arithmetic

Use native broadcasting and arithmetic with compatible cubes.

```julia
anom = fut .- ref
ratio = fut ./ ref
```

## Ensemble Helpers

```julia
stats = ensemble_stats(cube; dim="time")
stats2 = ensemble_fct(cube; dim="time")
```

`ensemble_fct` is an alias of `ensemble_stats`.
