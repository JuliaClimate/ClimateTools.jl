# Installation

## Stable Package Installation

Install ClimateTools from the general Julia package registry:

```julia
using Pkg
Pkg.add("ClimateTools")
```

ClimateTools targets Julia 1.10 and newer.

## Development Installation

To work on the package locally:

```julia
using Pkg
Pkg.develop(path="/path/to/ClimateTools")
```

## Core Dependencies

ClimateTools is built around:

- `YAXArrays.jl` for labeled gridded arrays and dataset loading
- `DimensionalData.jl` selectors for subsetting
- `Interpolations.jl`, `NearestNeighbors.jl`, and related tools for interpolation and regridding

Most workflows do not require a Python stack.

## Recommended Environment

For reproducible work, it is recommended to create a dedicated Julia environment:

```julia
using Pkg
Pkg.activate("climate-project")
Pkg.add(["ClimateTools", "YAXArrays", "Plots"])
```

## Optional Plotting Stack

ClimateTools itself does not require a specific plotting package. The documentation examples use generic Julia plotting tools rather than a package-specific plotting dependency.

Common choices include:

- `Plots.jl` for quick inspection plots
- `Makie.jl` for more customized figures
- GIS or Python tools outside Julia for publication-quality cartographic layouts

## Reading Datasets

Most users will also need `YAXArrays.jl` explicitly in their project when opening files:

```julia
using ClimateTools
using YAXArrays

cube = Cube(open_dataset("mydata.nc"))
```

## Building the Documentation Locally

If you want to build this documentation site from the repository root:

```julia
using Pkg
Pkg.activate("docs")
Pkg.instantiate()
include("docs/make.jl")
```

## Notes on Python and PyCall

Older ClimateTools workflows often relied on a Python-backed plotting stack. That is no longer required for the core package workflows documented here.

If your own visualization workflow uses `PyCall`, Cartopy, or other Python tools, manage that environment separately from the ClimateTools installation.

## Next Step

After installation, continue with the [Quick Start](quickstart.md) page.
