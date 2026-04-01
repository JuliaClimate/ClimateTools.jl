# ClimateTools.jl

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![CI](https://github.com/JuliaClimate/ClimateTools.jl/workflows/CI/badge.svg)](https://github.com/JuliaClimate/ClimateTools.jl/actions)
[![codecov](https://codecov.io/gh/JuliaClimate/ClimateTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaClimate/ClimateTools.jl)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaclimate.github.io/ClimateTools.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliaclimate.github.io/ClimateTools.jl/dev)
[![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/JuliaClimate/ClimateTools.jl)](https://github.com/JuliaClimate/ClimateTools.jl/releases)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/76293821.svg)](https://zenodo.org/badge/latestdoi/76293821)
[![chat](https://img.shields.io/badge/chat-on%20gitter-bc0067.svg)](https://gitter.im/ClimateTools-jl)

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
