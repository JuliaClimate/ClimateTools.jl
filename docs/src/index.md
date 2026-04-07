# ClimateTools.jl

ClimateTools.jl is a workflow-oriented climate-analysis package built on top of YAXArrays.jl. It is designed for users who need to move from raw gridded observations and simulations to usable climate scenarios, bias-corrected products, indicators, and summary diagnostics.

## What ClimateTools Covers

ClimateTools focuses on the common steps of applied climate analysis:

- Reading gridded CF-style datasets into `YAXArray` and `Cube` objects
- Selecting domains in time and space
- Regridding regular, curvilinear, and rotated-pole simulations
- Bias-correcting model simulations against observations or reanalyses
- Building climate scenarios from observational references and climate-model outputs
- Computing climate indicators, annual summaries, threshold counts, and spell statistics
- Performing selected diagnostics such as return-level estimation and time-series metrics

## Data Model

The package operates on YAXArrays-native objects.

```julia
using YAXArrays

cube = Cube(open_dataset("myfile.nc"))
```

High-level functions generally take a single `YAXArray` or a small group of aligned cubes and return another `YAXArray`, preserving axes where possible.

## Typical Workflow

A realistic ClimateTools workflow usually looks like this:

1. Open observational and simulation datasets.
2. Inspect dimensions, calendars, and coordinate conventions.
3. Subset the study period and spatial domain.
4. Interpolate or regrid datasets to a common grid.
5. Bias-correct model simulations.
6. Compute indicators and aggregations on the corrected data.
7. Validate and visualize the results.

```julia
using ClimateTools
using YAXArrays

obs = Cube(open_dataset("obs.nc"))
ref = Cube(open_dataset("ref.nc"))
fut = Cube(open_dataset("fut.nc"))

qq = qqmap(obs, ref, fut; method="additive", detrend=true)
txx = tx_max(qq)
```

## Start Here

If you are new to the package, read the pages in this order:

1. [Quick Start](quickstart.md)
2. [Data and Subsetting](datasets.md)
3. [Interpolation and Regridding](interpolation.md)
4. [Bias Correction](biascorrection.md)
5. [Building Climate Scenarios](scenarios.md)
6. [Indices and Aggregations](indices.md)

## Building Climate Scenarios

ClimateTools is especially useful for climate-scenario construction. In that context, a scenario usually means a climate-model projection that has been aligned to an observational reference through spatial interpolation and statistical post-processing.

The dedicated [Building Climate Scenarios](scenarios.md) guide walks through:

- choosing observational references
- preparing historical and future simulations
- regridding observations and simulations onto a common grid
- selecting a bias-correction method
- validating the resulting corrected fields
- computing derived indices for impact-oriented analysis

## Compute Model

Internally, ClimateTools prefers the YAXArrays `xmap` pattern for whole-dimension transforms and reductions such as time aggregations, period-wise indices, quantile summaries, and regridding over full grids.

Some workflows still use `mapCube`, mainly where multiple inputs share the same dimension name but not the same coordinate values. Bias-correction workflows based on observational, historical-model, and future-model series are the main example.

## Documentation Structure

- [Installation](installation.md): package setup and optional tooling
- [Quick Start](quickstart.md): first end-to-end result
- [Data and Subsetting](datasets.md): opening datasets, selecting ranges, polygon masking
- [Interpolation and Regridding](interpolation.md): regular and rotated-grid workflows
- [Bias Correction](biascorrection.md): method selection and usage guidance
- [Building Climate Scenarios](scenarios.md): full scenario workflow
- [Indices and Aggregations](indices.md): grouped catalog of indices and summaries
- [Examples](examples.md): worked examples by use case
- [Validation and Diagnostics](validation.md): how to assess outputs
- [Troubleshooting](troubleshooting.md): common pitfalls
- [API Overview](functions.md): grouped exported functions

## References and Attribution

Bias-correction and indicator workflows in ClimateTools are informed by the climate-services ecosystem, including xclim’s documentation style for workflow separation and the published literature behind each correction method. The references used in the docs are collected on the [References](references.md) page.
