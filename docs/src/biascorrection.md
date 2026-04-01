# Bias Correction

## Quantile Mapping

Use `qqmap` for day-of-year based quantile mapping on YAXArray/Cube data.

```julia
qq = qqmap(obs, ref, fut;
    method="additive",
    detrend=true,
    window=15,
    rankn=50,
    qmin=0.01,
    qmax=0.99)
```

## Bulk Variant

```julia
qq_bulk = qqmap_bulk(obs, ref, fut; method="multiplicative")
```

## Extreme-Value Tail Correction

Use `biascorrect_extremes` when you want multiplicative quantile mapping plus
an explicit GPD-based correction of upper-tail values.

```julia
dext = biascorrect_extremes(obs, ref, fut;
    detrend=false,
    P=0.95,
    runlength=2,
    frac=0.25,
    power=1.0)
```

The function accepts optional external extreme-value parameters with columns
`lat`, `lon`, `mu`, `sigma`, and `xi`:

```julia
using DataFrames

gevparams = DataFrame(
    lat=[45.0, 46.0],
    lon=[-73.0, -72.0],
    mu=[20.0, 21.5],
    sigma=[4.0, 4.5],
    xi=[0.10, 0.08],
)

dext = biascorrect_extremes(obs, ref, fut; gevparams=gevparams)
```
