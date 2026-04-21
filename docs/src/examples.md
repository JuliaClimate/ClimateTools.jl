# Examples

This page collects worked examples that connect multiple parts of the package instead of only isolated function calls.

## Example 1: Basic Bias Correction and Annual Maximum

```julia
using ClimateTools
using YAXArrays

obs = Cube(open_dataset("obs.nc"))
ref = Cube(open_dataset("ref.nc"))
fut = Cube(open_dataset("fut.nc"))

qq = qqmap(obs, ref, fut; method="additive", detrend=true)
txx = annualmax(qq)
```

Use this pattern when the goal is to derive a corrected future field and then summarize it annually.

## Example 2: Monthly Wet-Day Fraction

```julia
wet_fraction = wetdays_prop(pr; thresh=5.0, freq="MS")
```

This is a common first diagnostic when comparing present and future precipitation regimes.

## Example 3: Regrid Then Correct

```julia
regridder = Regridder(hist, obs; method="bilinear")
hist_rg = regrid(hist, regridder)
fut_rg = regrid(fut, regridder)

qq = qqmap(obs, hist_rg, fut_rg; method="additive")
```

This is the preferred order for many climate-scenario workflows: regrid first, then bias-correct.

## Example 4: Time Variability Correction

```julia
tvc_out = tvc(obs, hist, fut;
    scales=[365, 183, 92, 46, 23, 12, 6, 3, 2],
    eig_floor=1e-8)
```

Use this when correcting variability across time scales matters as much as correcting the marginal distribution.

## Example 5: Thermodynamic Post-Processing

```julia
vp = vaporpressure(huss, ps)
wb = wbgt(tas, vp)
```

This type of workflow is typical after you have already prepared a corrected scenario field.

## Example 6: Polygon-Based Regional Study

```julia
poly = ClimateTools.extractpoly("region.shp", n=1)
regional = spatialsubset(cube, poly)
regional_txx = tx_max(regional)
```

This pattern is useful for impact studies over watersheds, administrative regions, or user-defined polygons.

For rotated-pole or curvilinear outputs, use the parent dataset instead:

```julia
poly = ClimateTools.extractpoly("region.shp", n=1)
regional_ds = spatialsubset(ds, poly)
regional_txx = tx_max(regional_ds[:tasmax])
```

## Example 7: Faceted Geographic Comparison

```julia
fig = geomapfacet(fut;
    facetdim=:time,
    indices=1:4,
    ncols=2,
    colorbar=true,
    title="First four daily slices")
```

Use this when you want a quick multi-panel comparison over time with a shared projection and color range.

## Example 8: Ensemble Spread at One Grid Cell

```julia
fig = timeseriesplot(ensemble_cube;
    selectors=(longitude=10, latitude=12),
    mode=:mean_ribbon,
    spread=:minmax)
```

This is a compact diagnostic when comparing the temporal spread of multiple members at one location.

## Example 9: Grouped Ensemble Distributions

```julia
fig = statsplot(ensemble_cube;
    groupdim=:member,
    kind=:boxplot)
```

This is useful for checking whether one or two members dominate the overall ensemble spread.

## Example 10: Save and Reuse Extreme-Value Fits

```julia
using ClimateTools, YAXArrays, Serialization

cube = Cube(open_dataset("pr_daily.nc"))

gp_fit = gpfit_cube(cube;
    dim=:time,
    threshold_quantile=0.95,
    minimalvalue=0.0,
    return_thresholds=true,
    n_obs_per_block=365)

open("gpfit_cube.jls", "w") do io
    serialize(io, gp_fit)
end

gp_fit_loaded = open("gpfit_cube.jls", "r") do io
    deserialize(io)
end

rl = returnlevel_cube(gp_fit_loaded; rlevels=[2, 10, 25, 50, 100])
```

Use this pattern when fitting the gridpoint-wise tail model is expensive and you want to reuse the fitted GP models later without repeating the estimation step. The same reuse workflow applies to `gevfit_cube` after you have already reduced the source data to block maxima such as annual or seasonal maxima; it should not be applied directly to daily values. GP fits are the most complete example here because the saved bundle also carries the per-grid threshold cube needed for return-level calculations.

For a GEV workflow, reduce the source series to block maxima first and then fit across the resulting yearly maxima cube:

```julia
using ClimateTools, YAXArrays

cube = Cube(open_dataset("tasmax_daily.nc"))
annual_maxima = annualmax(cube)

gev_fit = gevfit_cube(annual_maxima; dim=:time)
gev_rl = returnlevel_cube(gev_fit; rlevels=[2, 10, 25, 50])
```

## Example 11: Monthly SPI and Water-Budget-First SPEI

```julia
using ClimateTools, Dates

spi3 = spi(pr;
    window=3,
    cal_start=DateTime(1981, 1, 1),
    cal_end=DateTime(2010, 12, 31))

water_budget = pr .- pet

spei6 = spei(water_budget;
    window=6,
    cal_start=DateTime(1981, 1, 1),
    cal_end=DateTime(2010, 12, 31))
```

This first pass uses a gamma-based monthly core. `spi` applies zero inflation for zero-precipitation months, while `spei` assumes the supplied water-budget series is already prepared and positive for the period being fit.

## Where to Go Next

- [Quick Start](quickstart.md) for the minimal onboarding path
- [Building Climate Scenarios](scenarios.md) for the full observational-plus-simulation workflow
- [Validation and Diagnostics](validation.md) for how to assess outputs
