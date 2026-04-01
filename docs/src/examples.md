# Examples

## End-to-end Workflow

```julia
using ClimateTools
using YAXArrays

obs = Cube(open_dataset("obs.nc"))
ref = Cube(open_dataset("ref.nc"))
fut = Cube(open_dataset("fut.nc"))

qq = qqmap(obs, ref, fut; method="additive", detrend=true)
ann = annualmax(qq)

# Optional regrid
# ann_rg = regrid_cube(ann, target_cube)
```

## Daily Aggregation

```julia
daily_mean = daymean(fut)
daily_sum = daysum(fut)
```

## Thermodynamic Helpers

```julia
vp = vaporpressure(huss, ps)
wb = wbgt(tas, vp)
```
