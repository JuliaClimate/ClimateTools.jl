# Visualization

ClimateTools now includes a small plotting API for the most common climate-analysis views:

- `geomap` for geographic fields
- `geomapfacet` for multi-panel geographic comparisons over time or ensemble dimensions
- `timeseriesplot` for point or reduced time series
- `statsplot` for quick distribution and comparison plots

These methods are available when `GeoMakie.jl` and a Makie backend are loaded.

## Setup

```julia
using ClimateTools
using GeoMakie
using CairoMakie

CairoMakie.activate!()
```

## Geographic Maps

For regular longitude-latitude grids, `geomap` can infer the spatial axes automatically.

```julia
fig = geomap(cube; dim=:time, index=1)
```

For a more publication-oriented layout, `geomap` now exposes explicit map controls for projection, framing, and colorbars:

```julia
fig = geomap(cube;
    dim=:time,
    index=1,
    colorbar=true,
    colorbar_label="Temperature (K)",
    colorbar_position=:right,
    coastline_width=1.5,
    dest="+proj=eqearth +lon_0=0",
    frame=false)
```

The most useful map keywords are:

- `source`, `dest`, `lon_0` for explicit projection control
- `limits`, `fit_limits`, `limit_padding` for controlling the displayed longitude-latitude extent
- `colorbar`, `colorbar_label`, `colorbar_position` for layout and legend-like scaling
- `coastline_color`, `coastline_width` for cartographic emphasis
- `frame` to suppress axis decorations in cleaner publication layouts
- `axis_kwargs`, `surface_kwargs`, `colorbar_kwargs` for lower-level Makie customization

By default, `geomap` and `geomapfacet` fit the visible map extent to the minimum and maximum longitude and latitude of the plotted data. This works for both ordinary lon-lat grids and rotated-pole inputs once they have been transformed back to geographic coordinates.

If you want to override the default view window explicitly:

```julia
fig = geomap(cube;
    dim=:time,
    index=1,
    limits=((-75.0, -55.0), (43.0, 63.0)))
```

If you want GeoMakie to manage the limits instead, disable the automatic fit:

```julia
fig = geomap(cube; dim=:time, index=1, fit_limits=false)
```

If your cube still has additional non-spatial dimensions, use `selectors` to reduce them first or specify which dimension should be sliced:

```julia
fig = geomap(cube; selectors=(member=1,), dim=:time, index=10)
```

For rotated-pole or curvilinear datasets, pass the parent dataset and variable name so ClimateTools can recover the geographic coordinates from `grid_mapping`, `lon`, and `lat` metadata:

```julia
ds = open_dataset("rotated_model.nc")
fig = geomap(ds, :tas; dim=:time, index=1)
```

The default projection is selected from the coordinate extent. Global fields use an equal-earth style projection, while regional fields default to geographic longitude-latitude coordinates.

## Faceted Maps

Use `geomapfacet` when you want several aligned panels with a shared projection and shared color scale.

```julia
fig = geomapfacet(cube;
    facetdim=:time,
    ncols=3,
    colorbar=true,
    title="Daily fields")
```

If you want all panels to use one common geographic window, even when the facet payloads do not share the same spatial extent, enable `shared_spatial_limits`:

```julia
fig = geomapfacet(cube;
    facetdim=:time,
    ncols=3,
    shared_spatial_limits=true)
```

This is separate from `limits=...`:

- `shared_spatial_limits=true` computes one unioned longitude-latitude window from all facet payloads
- `limits=((lonmin, lonmax), (latmin, latmax))` enforces an explicit user-defined window

Projection sharing is now explicit as well:

- `shared_projection=true` keeps one common destination projection across all panels
- `shared_projection=false` allows `dest` to vary per panel
- `shared_dest` is accepted as an alias for the same boolean behavior

With a shared projection, pass one `dest` string for all panels:

```julia
fig = geomapfacet(cube;
    facetdim=:time,
    indices=1:2,
    shared_projection=true,
    dest="+proj=eqearth +lon_0=0")
```

If you want panel-specific projections, disable projection sharing and pass either a collection or a function:

```julia
fig = geomapfacet(cube;
    facetdim=:time,
    indices=1:2,
    shared_projection=false,
    dest=["+proj=eqearth +lon_0=0", "+proj=longlat +datum=WGS84"])
```

For ensemble data, facet over the ensemble dimension and fix the time slice with `selectors`:

```julia
fig = geomapfacet(cube;
    facetdim=:member,
    selectors=(time=1,),
    ncols=2)
```

`geomapfacet` is especially useful when comparing:

- sequential dates from the same simulation
- ensemble members at a fixed date
- raw and corrected outputs after regridding or bias correction

The default behavior uses a shared color range across panels. This keeps visual comparisons honest, especially for multi-member diagnostics.

## Time-Series Inspection

Extract a single grid-cell series with selectors:

```julia
fig = timeseriesplot(cube; selectors=(longitude=10, latitude=12))
```

Or collapse the remaining spatial dimensions with a reducer such as `mean`:

```julia
fig = timeseriesplot(cube; reducer=mean, ylabel="Regional mean")
```

When a cube contains an ensemble-like dimension such as `member`, `scenario`, or `run`, `timeseriesplot` can render all members directly:

```julia
fig = timeseriesplot(cube;
    selectors=(longitude=10, latitude=12),
    mode=:lines)
```

Or collapse the ensemble into a central line plus spread ribbon:

```julia
fig = timeseriesplot(cube;
    selectors=(longitude=10, latitude=12),
    mode=:mean_ribbon,
    spread=:minmax)
```

If you already summarized a cube with `ensemble_stats`, you can visualize the resulting `:stats` dimension directly:

```julia
stats = ensemble_stats(cube; dim="member")
fig = timeseriesplot(stats; selectors=(longitude=10, latitude=12), mode=:stats)
```

You can also compare several already extracted series directly:

```julia
fig = timeseriesplot((obs=obs_series, sim=sim_series, corrected=bc_series))
```

## Statistical Summaries

`statsplot` is intended for fast diagnostics rather than publication styling.

```julia
fig = statsplot((raw=raw_values, corrected=corrected_values); kind=:hist)
fig = statsplot((raw=raw_values, corrected=corrected_values); kind=:boxplot)
fig = statsplot((obs=obs_values, sim=sim_values); kind=:scatter)
```

For a single vector or cube slice, histogram and boxplot views are supported directly.

For ensemble cubes, you can group distributions by member or scenario:

```julia
fig = statsplot(cube; groupdim=:member, kind=:boxplot)
```

## When to Use These Helpers

These functions are meant for:

- workflow diagnostics after regridding or bias correction
- comparing raw and corrected simulations
- checking temporal behavior at a grid cell or after regional averaging
- quickly inspecting rotated-grid outputs without first hand-building coordinate transforms

For highly customized figure layouts, annotations, or publication graphics, it is still reasonable to extract arrays from ClimateTools outputs and build Makie figures manually.

## What to Plot During Scenario Workflows

Useful visual checks include:

- the reference observation field
- the raw model field on the same date or period
- the regridded model field
- the corrected field
- the difference between corrected and raw simulations

For bias-correction work, plotting a few representative time series is often as informative as map slices.

## Related Pages

- [Examples](examples.md)
- [Validation and Diagnostics](validation.md)
- [Building Climate Scenarios](scenarios.md)
