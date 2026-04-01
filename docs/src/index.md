# ClimateTools.jl

ClimateTools.jl provides climate-analysis functions on top of YAXArrays.jl.

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
