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

ClimateTools.jl is a climate-analysis toolkit built on top of YAXArrays.jl for gridded observations, reanalyses, and climate-model simulations.

The package focuses on the typical climate-services workflow:

1. Read CF-style gridded datasets as `YAXArray` or `Cube` objects.
2. Subset them in time and space.
3. Regrid simulations and observations onto a common domain.
4. Bias-correct simulations with quantile mapping, extreme-value correction, or time-variability correction.
5. Compute indicators, aggregations, and diagnostics.

ClimateTools switched to YAXArrays `xmap` and `mapCube` pattern for whole-dimension transforms and reductions instead of the custom `ClimGrid`.

## Status

- Active development
- Julia 1.10+
- Data model: `YAXArray` / `Cube`
- YAXArrays compatibility: `0.7`

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

## Quick Start

```julia
using ClimateTools
using YAXArrays

obs = Cube(open_dataset("obs.nc"))
ref = Cube(open_dataset("ref.nc"))
fut = Cube(open_dataset("fut.nc"))

qq = qqmap(obs, ref, fut; method="additive", detrend=true)
ann = annualmax(qq)
```

That small example already illustrates the main ClimateTools workflow: open data, bias-correct, then compute a climate summary.

For quantile mapping, `qqmap` is the seasonally varying method: it removes leap days, groups samples by day of year using a moving `+/- window`, and applies an additive or multiplicative correction that changes through the annual cycle. Use `qqmap_bulk` only when one full-sample correction is acceptable and seasonal variation in the bias is not important.

## Main Capabilities

- Dataset opening and subsetting with YAXArrays and DimensionalData selectors
- Regridding on regular, curvilinear, and rotated-pole grids
- Bias correction with `qqmap`, `qqmap_bulk`, `biascorrect_extremes`, and `tvc`
- Aggregations such as `daymean`, `daysum`, `annualmax`, `annualmean`, and `annualsum`
- A large catalog of temperature, precipitation, threshold, and spell-duration indices
- Ensemble summaries and time-series diagnostics
- Thermodynamic helper functions and return-level estimation

## Documentation Map

The stable documentation is organized around workflows first and reference second.

- Quick start: open data, inspect dimensions, compute a first result
- Data and subsetting: time ranges, polygons, coordinate conventions
- Interpolation and regridding: regular, rotated, and curvilinear workflows
- Bias correction: method selection, assumptions, and validation
- Building climate scenarios: observations plus simulations, interpolation, correction, diagnostics, and derived indicators
- Indices and aggregations: grouped by climate-analysis task
- API overview: exported functions by category

## Climate Scenarios

One of the primary use cases for ClimateTools is climate-scenario construction. A typical scenario workflow is:

1. Choose a reference observation or reanalysis dataset.
2. Open one or more model simulations for a historical and future period.
3. Align calendars, time windows, and grids.
4. Interpolate observations and simulations to a common grid if necessary.
5. Bias-correct the model series using a method suited to the variable and objective.
6. Validate the corrected fields and compute indicators, maps, or summary tables.

The documentation includes a dedicated guide for this workflow.

## Documentation

See the built documentation at:

- Stable: [https://juliaclimate.github.io/ClimateTools.jl/stable](https://juliaclimate.github.io/ClimateTools.jl/stable)
- Dev: [https://juliaclimate.github.io/ClimateTools.jl/dev](https://juliaclimate.github.io/ClimateTools.jl/dev)
