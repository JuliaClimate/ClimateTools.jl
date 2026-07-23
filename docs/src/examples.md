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

For a stationary precipitation GEV workflow, reduce daily precipitation amounts to annual maximum one-day precipitation before fitting across the block axis:

```julia
using ClimateTools, YAXArrays

pr = Cube(open_dataset("pr_daily.nc"))
annual_max_precipitation = annualmax(pr)

gev_fit = gevfit_cube(annual_max_precipitation; dim=:time)
gev_rl = returnlevel_cube(gev_fit; rlevels=[2, 10, 25, 50])
```

For a non-stationary GEV model, a greenhouse-gas concentration can be supplied along the block axis. This example uses the [NOAA Global Monitoring Laboratory global annual mean atmospheric CO2 record](https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_annmean_gl.csv), expressed in ppm. A snapshot covering 1979–2025 is bundled with the package test data. The concentration records are matched by year rather than by row position, and centering the covariate improves numerical conditioning:

```julia
using ClimateTools, YAXArrays, Dates, Statistics, DelimitedFiles

pr = Cube(open_dataset("pr_daily.nc"))
annual_max_precipitation = annualmax(pr)

co2_file = joinpath(pkgdir(ClimateTools), "test", "data", "co2_annmean_gl.csv")
co2_data, _ = readdlm(co2_file, ',', Float64;
    header=true,
    comments=true,
    comment_char='#')

co2_by_year = Dict(Int(row[1]) => row[2] for row in eachrow(co2_data))
block_years = year.(lookup(annual_max_precipitation, :time))
co2_ppm = [get(co2_by_year, block_year, missing) for block_year in block_years]
co2_reference = mean(skipmissing(co2_ppm))
centered_co2_ppm = [
    ismissing(value) ? missing : value - co2_reference
    for value in co2_ppm
]

gev_fit = gevfit_cube(annual_max_precipitation;
    dim=:time,
    covariates=(co2_ppm=centered_co2_ppm,),
    locationcovid=[:co2_ppm])

model = first(skipmissing(Array(gev_fit.models)))
coefficients = model.θ̂
effective_rl = returnlevel_cube(gev_fit; rlevels=[10, 25, 50, 100])
```

Here `coefficients` contains the location intercept and the change in location per ppm of atmospheric CO2, followed by the stationary log-scale and shape coefficients. `Extremes.params(model)` can be used when the location, scale, and shape values evaluated at every retained covariate row are needed instead. The effective return-level cube restores the fitted `time` axis, so its dimensions are the remaining input dimensions followed by `time` and `rlevels`. Years outside the NOAA record, along with rows containing `missing`, `NaN`, or infinite precipitation or covariate values, remain `missing` in that output. Each grid cell must have at least `min_points` complete rows and at least as many complete rows as fitted regression coefficients.

Atmospheric CO2 is one greenhouse-gas covariate, not a combined concentration of all greenhouse gases. For future climate simulations, replace the observational NOAA series with concentrations from the same emissions scenario as the precipitation simulation. The precipitation cube should contain daily amounts or totals in consistent units; convert rates or fluxes before extracting annual maxima.

Covariates may also be YAXArrays. A covariate cube must include the fitted dimension and may include any subset of the remaining dimensions for spatially or ensemble-varying regressions. Coordinates on every shared axis must match the block-maxima cube exactly. The `locationcovid`, `logscalecovid`, and `shapecovid` keywords follow the corresponding Extremes.jl non-stationary GEV interface.

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
