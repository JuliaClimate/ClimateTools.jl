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

For xclim-style ensemble post-processing, ClimateTools now also provides richer summary and reduction helpers around a `realization` axis.

```julia
ens_stats = ensemble_mean_std_max_min(cube; realization_dim="realization")
ens_perc = ensemble_percentiles(cube;
    realization_dim="realization",
    values=[10, 50, 90])
```

`ensemble_mean_std_max_min` returns a dataset of mean, standard deviation, maximum, and minimum cubes.

`ensemble_percentiles` returns split percentile variables such as `p10`, `p50`, and `p90` by default. Use `split=false` to keep the percentiles on a dedicated `percentiles` dimension.

Weighted summaries are supported in both helpers through the `weights` keyword.

These functions are useful when:

- working with multi-member simulations
- summarizing a stack of realizations
- deriving ensemble-level diagnostics after bias correction

The reduction layer includes deterministic preprocessing and subset selection:

```julia
criteria = make_criteria(cube; realization_dim="realization")
selected = kkz_reduce_ensemble(criteria, 5)
```

`make_criteria` flattens all non-realization dimensions into a single `criteria` axis.

`kkz_reduce_ensemble` applies the KKZ algorithm to select a representative subset of members without introducing stochastic clustering dependencies.

Robustness diagnostics for change signals are also available:

```julia
fractions = robustness_fractions(fut, ref; test="threshold", rel_thresh=0.02)
classes = robustness_categories(fractions)
R = robustness_coefficient(fut, ref)
```

`robustness_fractions` summarizes how much of the ensemble changes, in which direction, and how strongly members agree.

`robustness_categories` converts those fractions into IPCC-style robustness classes.

`robustness_coefficient` implements the Knutti-Sedlacek robustness coefficient.

This first ensemble port focuses on xclim-compatible summaries, criteria flattening, deterministic KKZ reduction, and robustness metrics. K-means reduction and partitioning methods are still separate follow-up work.

Those summaries now connect directly to the plotting layer. Two common patterns are:

```julia
fig = timeseriesplot(cube;
    selectors=(longitude=10, latitude=12),
    mode=:mean_ribbon)

stats = ensemble_stats(cube; dim="member")
fig2 = timeseriesplot(stats;
    selectors=(longitude=10, latitude=12),
    mode=:stats)

xclim_stats = ensemble_mean_std_max_min(cube; realization_dim="member")
fig3 = timeseriesplot(xclim_stats, :mean;
    selectors=(longitude=10, latitude=12))

xclim_percentiles = ensemble_percentiles(cube;
    realization_dim="member",
    values=[10, 50, 90],
    split=false)
fig4 = timeseriesplot(xclim_percentiles;
    selectors=(longitude=10, latitude=12))
```

`timeseriesplot(dataset, :varname)` and `statsplot(dataset, :varname)` now work directly on dataset outputs such as `ensemble_mean_std_max_min`.

When `ensemble_percentiles(...; split=false)` is used, `timeseriesplot` automatically recognizes the `percentiles` dimension and draws one line per requested percentile.

For geographic comparison of members, use `geomapfacet(cube; facetdim=:member, selectors=(time=1,))`.

For robustness diagnostics, the plotting layer also includes a categorical map helper:

```julia
fractions = robustness_fractions(fut, ref; test="threshold", rel_thresh=0.02)
fig5 = robustnessmap(fractions)

classes = robustness_categories(fractions)
fig6 = robustnessmap(classes)
```

`robustnessmap` accepts either the output of `robustness_fractions` or a precomputed categorical result from `robustness_categories`, and renders a discrete category map with labels derived from the robustness metadata.

## Where This Fits in the Workflow

Arithmetic and ensemble summaries are often the final layer after you have:

1. opened and aligned the data
2. regridded the datasets if needed
3. bias-corrected the simulation fields

At that stage, cube arithmetic becomes a compact way to compare scenarios.
