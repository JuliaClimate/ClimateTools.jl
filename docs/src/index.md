# ClimateTools.jl

ClimateTools.jl provides climate-analysis functions on top of YAXArrays.jl.

ClimateTools now uses the YAXArrays `xmap` pattern for most whole-dimension reductions and transforms while preserving the same public YAXArray-oriented API.

## Data Model

All high-level APIs operate on YAXArray/Cube objects.

```julia
using YAXArrays
cube = Cube(open_dataset("myfile.nc"))
```

## Typical Pipeline

1. Read data as YAXArray/Cube.
2. Apply preprocessing or bias correction.
3. Compute daily/yearly summaries and indices.
4. Regrid or visualize results.

```julia
using ClimateTools
qq = qqmap(obs, ref, fut; method="additive")
ann = annualmax(qq)
```

## Compute Pattern

If you are extending ClimateTools internals, prefer `xmap` for operations that reduce or replace complete dimensions such as time aggregations, period reductions, quantile summaries, or regridding over full source grids.

`mapCube` still appears in a few specialized implementations where the current YAXArrays `xmap` rules are more restrictive. The main case is a multi-input workflow where arrays share a dimension name but do not share the same coordinate values, such as the `obs`/`ref`/`fut` time axes used by quantile mapping.

## Included Documentation

- Installation
- Datasets and subsetting
- Bias correction
- Indices and aggregations
- Interpolation and regridding
- Examples
- Function overview (full exported API list)

## Additional Capabilities

- Extreme-value bias correction with `biascorrect_extremes`
- Return-level estimation with `rlevels_cube`
- Time-series diagnostics with `autocorrelation` and `hurst`
