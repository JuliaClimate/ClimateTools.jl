# ClimateTools.jl

ClimateTools.jl is a climate-analysis toolkit built on top of YAXArrays.jl.

## Status

- Active development
- Julia 1.10+
- Data model: YAXArray/Cube

## Installation

```julia
using Pkg
Pkg.add("ClimateTools")
```

For local development:

```julia
using Pkg
Pkg.develop(path="/path/to/ClimateTools")
```

## Core Workflow (YAXArrays-native)

```julia
using ClimateTools
using YAXArrays

obs = Cube(open_dataset("obs.nc"))
ref = Cube(open_dataset("ref.nc"))
fut = Cube(open_dataset("fut.nc"))

qq = qqmap(obs, ref, fut; method="additive", detrend=true)
ann = annualmax(qq)
```

## Main Features

- Quantile mapping bias correction: `qqmap`, `qqmap_bulk`
- Daily and yearly aggregation: `daily_fct`, `daymean`, `daysum`, `yearly_resample`, `monthly_resample`
- Annual indices and threshold counts: `annualmax`, `annualmin`, `annualmean`, `annualsum`, `prcp1`, `frostdays`, `summerdays`, `icingdays`, `tropicalnights`, `customthresover`, `customthresunder`
- Thermodynamic helpers: `vaporpressure`, `approx_surfacepressure`, `wbgt`, `diurnaltemperature`, `meantemperature`
- Regridding utilities: `regrid_cube`, `regrid_curvilinear_to_regular`, `regrid_rotated_curvilinear_to_regular`
- Ensemble summaries: `ensemble_stats`, `ensemble_fct`

## Documentation

See `docs/src/` for full YAXArrays-based usage examples.
