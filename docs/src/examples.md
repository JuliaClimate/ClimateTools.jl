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

## Where to Go Next

- [Quick Start](quickstart.md) for the minimal onboarding path
- [Building Climate Scenarios](scenarios.md) for the full observational-plus-simulation workflow
- [Validation and Diagnostics](validation.md) for how to assess outputs
