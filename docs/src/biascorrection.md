# Bias Correction

ClimateTools provides several bias-correction methods aimed at different goals. Some methods mainly correct distributions, while others target time-scale variability or upper-tail behavior.

Before choosing a method, make sure that:

- the observational reference and the model simulation cover a comparable calibration period
- the datasets have compatible units
- the grids are compatible, or have already been regridded to a common target

## Which Method Should You Use?

Use this quick guide:

- `qqmap`: default day-of-year quantile mapping for many temperature and precipitation workflows
- `qqmap_bulk`: bulk quantile mapping without day-of-year grouping
- `biascorrect_extremes`: upper-tail-aware correction when extremes matter
- `tvc`: correction of temporal covariance across time scales while preserving event order

## Quantile Mapping with `qqmap`

`qqmap` performs day-of-year-based quantile mapping.

```julia
qq = qqmap(obs, ref, fut;
    method="additive",
    detrend=true,
    window=15,
    rankn=50,
    qmin=0.01,
    qmax=0.99)
```

Typical choices:

- `method="additive"` for temperature-like variables
- `method="multiplicative"` for precipitation-like variables

Important parameters:

- `window`: day-of-year neighborhood used for seasonal grouping
- `rankn`: number of quantiles used in the mapping
- `detrend`: whether to remove and reapply a polynomial trend during correction

Use `qqmap` when seasonal behavior matters and you want the correction to vary through the annual cycle.

### How `qqmap` Works

`qqmap` compares three series:

- `obs`: the observed reference
- `ref`: the model series over the calibration period
- `fut`: the model series to correct

The correction is computed grid point by grid point. ClimateTools first removes February 29 from all three inputs so that every year is treated as a 365-day calendar. This avoids mixing leap-day values into the seasonal grouping.

For each Julian day, ClimateTools:

1. selects the `fut` values that occur on that exact day of year
2. gathers `obs` and `ref` values from a moving window of `+/- window` days around that day across all years
3. estimates empirical quantiles for the `obs` and `ref` window samples
4. builds an interpolation from the `ref` quantiles to a correction term
5. applies that correction to the matching `fut` values for that day

This is why `qqmap` is seasonal without using explicit monthly bins. The correction evolves smoothly through the annual cycle because each day of year gets its own transfer function, estimated from a local day-of-year neighborhood. Near New Year, the moving window wraps around the year boundary, so late-December and early-January samples are grouped together rather than truncated.

### Additive and Multiplicative Branches

ClimateTools implements two correction formulas:

- additive: build the quantile-wise shift `obsP - refP`, then apply `x_corr = x_fut + f(x_fut)`
- multiplicative: build the quantile-wise scale factor `obsP / refP`, then apply `x_corr = x_fut * f(x_fut)`

In practice:

- use additive correction for temperature-like variables where a shift is physically meaningful
- use multiplicative correction for precipitation-like variables where relative scaling is more appropriate

Interpolation is built on the estimated `ref` quantiles, and the default extrapolation is flat at the lower and upper ends. That means values beyond the fitted quantile range receive the edge correction rather than an unbounded extrapolated one.

### Detrending Behavior

If `detrend=true`, ClimateTools fits a polynomial trend to each input series as a function of time index, removes that fitted trend before quantile mapping, and then restores the fitted trend of `fut` after correction. The intention is to correct the distribution of anomalies while preserving the large-scale trend already present in the climate simulation.

### When to Use `qqmap` Versus `qqmap_bulk`

Choose `qqmap` when biases depend on the season, for example when winter and summer model errors differ or when wet-season and dry-season precipitation biases should not share the same correction.

Choose `qqmap_bulk` only when a single distributional correction over the full sample is acceptable. `qqmap_bulk` uses the full time series at once and does not preserve seasonal variation in the correction. It is simpler, but that simplicity can smear seasonally varying biases into one annual mapping. *It is useful for weather prediction bias correction*

### Practical Caveats

- If too many values in the local `obs`, `ref`, or `fut` samples are `NaN`, the correction is skipped for that location and time, returning either `NaN` or the original future values depending on `keep_original` keyword (default to `false`).
- Multiplicative correction can still be fragile near zero or with sparse wet-day samples.
- The moving-window approach improves seasonal realism, but very short calibration series can still produce unstable quantiles for some days of year.

## Bulk Quantile Mapping with `qqmap_bulk`

```julia
qq_bulk = qqmap_bulk(obs, ref, fut; method="multiplicative")
```

This method uses the full sample rather than a moving day-of-year window. It is simpler, but it does not adapt the correction seasonally.

## Extreme-Value Tail Correction with `biascorrect_extremes`

Use `biascorrect_extremes` when the upper tail needs special treatment beyond standard quantile mapping. **Tested on precipitation fields**

This function is the ClimateTools implementation of the mixed QQM-GPD strategy described by Roy et al. (2023), which combines non-parametric quantile mapping for the bulk of the precipitation distribution with an extreme-value tail correction for high precipitation. See [References](references.md).

```julia
dext = biascorrect_extremes(obs, ref, fut;
    detrend=false,
    P=0.95,
    runlength=2,
    frac=0.25,
    power=1.0)
```

The function combines multiplicative quantile mapping with an explicit generalized Pareto distribution correction of extreme values.

### How `biascorrect_extremes` Is Implemented

`biascorrect_extremes` is a two-stage method.

1. ClimateTools first computes a base field using `qqmap(obs, ref, fut; method="multiplicative")`.
2. It then revisits only the upper-tail events and selectively replaces those values using a generalized Pareto distribution (GPD) mapping.

The output is initialized from the multiplicative `qqmap` result. If an extreme-tail correction cannot be estimated at a grid point, the method keeps the base `qqmap` values rather than failing the whole field.

In the paper, this mixed method is described as a QQM-GPD approach: QQM handles the bulk distribution and a parametric generalized Pareto tail model handles the right tail, with a smooth transition between the two.

### Mathematical Formulation

In the notation of Roy et al. (2023), the corrected precipitation distribution is written as a semi-parametric mixture between a bulk distribution $H_X$ and a right-tail GPD term $G_{(X-\ell \mid X > \ell)}$:

$$
F_X(x; h, k, r, n)
=
\left(1 - w(x; k)\right) H_X(x; h)
+ w(x; k) G_{(X-\ell \mid X > \ell)}(x - \ell; r, n).
$$

Here:

- $H_X$ is the bulk distribution estimated with QQM
- $G$ is the generalized Pareto tail model above threshold $\ell$
- $w(x; k)$ is a transition weight increasing from $0$ in the bulk to $1$ in the far right tail

The paper defines a lower threshold $\ell$ and an upper transition threshold $t$, with a linear transition between them:

$$
w(x; \ell, t) =
\begin{cases}
0, & x \le \ell, \\
\dfrac{H_X(x; h) - H_X(\ell; h)}{H_X(t; h) - H_X(\ell; h)}, & \ell < x < t, \\
1, & x \ge t.
\end{cases}
$$

For precipitation bias correction, the QQM mapping is then written as

$$
F_{Y_F}(x) = F_{Y_C}\!\left(F_{X_C}^{-1}\!\left(F_{X_F}(x)\right)\right),
$$

where $Y_C$ is observed current precipitation, $X_C$ is simulated current precipitation, and $X_F$ is simulated future precipitation. In the mixed QQM-GPD method, this QQM relation is used for the bulk, while the right tail is corrected in GPD probability space before being blended back into the full distribution.

In the current ClimateTools implementation, the paper symbols map approximately as follows: $H_X$ corresponds to the multiplicative `qqmap` baseline used as the bulk correction; $G$ corresponds to the tail models returned by `_fit_gpd_distribution` or `estimate_gpd`; the lower threshold $\ell$ corresponds to the threshold computed by `_threshold_from_vectors` and controlled by `P`; the upper transition threshold $t$ is not stored explicitly but is represented operationally by `_extremes_weight` through the keywords `frac` and `power`; and the transition weight $w$ is the blend produced by `_extremes_weight`, which determines how strongly the final corrected extremes depart from the `qqmap` baseline toward the GPD-based tail correction.

### Core Steps

For each grid point, the implementation does the following:

1. infer the time dimensions of `obs`, `ref`, and `fut`
2. drop February 29 from all three inputs
3. compute the multiplicative `qqmap` base correction
4. flatten the spatial dimensions into independent time series
5. skip locations where the `NaN` fraction in `obs`, `ref`, or `fut` exceeds `thresnan`
6. compute an extreme threshold from the mean of the `P` quantiles of `obs` and `ref`, using only finite values greater than or equal to `1.0`
7. fit or retrieve an observed-tail GPD
8. fit a future-tail GPD on clustered threshold exceedances
9. map future exceedance probabilities from the future GPD into the observed-tail GPD
10. blend the tail-corrected values back into the multiplicative `qqmap` baseline

### Threshold And Clustering Logic

The threshold is not estimated from `fut`. ClimateTools computes it from the calibration data using:

- the `P` quantile of finite `obs` values above `1.0`
- the `P` quantile of finite `ref` values above `1.0`
- the mean of those two quantiles

That threshold is then used to detect extreme clusters with `Extremes.getcluster`. The `runlength` keyword controls the cluster definition. For each cluster, ClimateTools keeps the cluster maximum and fits a GPD to the exceedances above the threshold.

This makes the tail correction event-oriented rather than pointwise over every time step above the threshold.

This is consistent with the paper’s motivation: heavy precipitation is treated as clustered meteorological events rather than as independent daily exceedances. The paper uses a high-tail threshold at the 95th percentile and a 1.0 mm event-separation rule for declustering. In ClimateTools, the analogous controls are `P=0.95` by default and the clustering logic implemented through `Extremes.getcluster`.

### Internal Versus External Tail Models

There are two branches for the observed-tail distribution.

- If `gevparams` is empty, ClimateTools estimates the observed-tail GPD directly from clustered exceedances in `obs`.
- If `gevparams` is provided, ClimateTools finds the nearest row in `gevparams` using the local latitude and longitude, then converts the provided GEV parameters `mu`, `sigma`, and `xi` into a GPD-like tail model at the chosen threshold.

If `gevparams` is supplied, the input cube must expose latitude and longitude dimensions that ClimateTools can infer. Otherwise the function errors because it cannot match grid points to parameter rows.

The external-parameter branch directly reflects the paper’s intended workflow for hydrological applications: use externally estimated GEV parameters from a more robust extreme-value dataset, then convert them to a GPD tail model at the chosen threshold before bias correcting the simulated extremes.

### Long-Series Moving Windows

If the future period is much longer than the reference period, ClimateTools switches to a moving-window tail correction. Concretely, this happens when the year span of `fut` is greater than `1.5` times the year span of `ref`.

In that case, ClimateTools:

- slices `fut` into overlapping windows whose length is based on `ref`
- fits the future-tail GPD separately in each window
- applies the tail correction only in the central portion of each window

This avoids forcing one tail model estimated over a very long future series to represent all periods equally.

### Probability Mapping And Blending

Once both future and observed-tail GPDs are available, ClimateTools:

1. computes the CDF of the future exceedances under the future-tail GPD
2. clamps those probabilities to the interval `[1e-6, 1 - 1e-6]`
3. evaluates the corresponding quantiles under the observed-tail GPD
4. adds the threshold back to obtain corrected extreme values

The corrected extremes are not inserted as a hard replacement. ClimateTools computes a transition weight from the exceedance magnitude using `frac` and `power`, then blends the tail-corrected values with the multiplicative `qqmap` baseline.

That blend is important: lower exceedances remain closer to the base `qqmap` field, while the strongest exceedances move more fully toward the tail-corrected values.

In Roy et al. (2023), the published QQM-GPD formulation uses a linear transition between a lower threshold $l$ and an upper threshold $t = l + \frac{1}{4}(\max(X) - l)$. The current ClimateTools implementation is aligned with that idea but implements the blend through `frac` and `power` in `_extremes_weight`. With the defaults `P=0.95`, `frac=0.25`, and `power=1.0`, the code follows the same practical intent: leave non-extreme values on the QQM baseline and transition progressively toward the parametric tail correction as exceedance magnitude increases.

### Practical Interpretation Of The Keywords

- `P`: controls the threshold used to define the extreme tail
- `runlength`: controls the temporal clustering of exceedances
- `gevparams`: switches from locally fitted observed tails to externally supplied tail parameters
- `frac` and `power`: control how aggressively the tail correction replaces the base `qqmap` values
- `detrend`: propagates into the base multiplicative `qqmap` stage before the tail adjustment is applied

### Caveats

- The base correction is always multiplicative `qqmap`, so this method is mainly aligned with precipitation-like or strictly positive upper-tail variables.
- If no stable GPD fit is available at a grid point, the method silently falls back to the multiplicative `qqmap` baseline at that location.
- The threshold calculation ignores values below `1.0`, which is a meaningful modeling choice for precipitation-like data but may be inappropriate for all variables.
- The internal GPD fitting path uses `try/catch` and can skip problematic locations without raising an error.
- The current implementation is paper-aligned but not a literal transcription of the mixture-CDF equations. It reproduces the same workflow idea, but the transition is implemented through value-based blending on corrected extremes rather than by explicitly evaluating the semi-parametric mixture distribution written in the article.

The function also accepts optional external extreme-value parameters with columns `lat`, `lon`, `mu`, `sigma`, and `xi`:

```julia
using DataFrames

gevparams = DataFrame(
    lat=[45.0, 46.0],
    lon=[-73.0, -72.0],
    mu=[20.0, 21.5],
    sigma=[4.0, 4.5],
    xi=[0.10, 0.08],
)

dext = biascorrect_extremes(obs, ref, fut; gevparams=gevparams)
```

This method is most appropriate when impact-relevant extremes are central to the study.

## Time Variability Correction with `tvc`

`tvc` applies the Time Variability Correction method of Shao et al. (2024), DOI 10.1029/2023MS003640.

The method is designed to correct covariance across and between multiple time scales while preserving the event sequence of the validation series.

```julia
corrected = tvc(obs, raw_train, raw_val;
    scales=[365, 183, 92, 46, 23, 12, 6, 3, 2],
    eig_floor=1e-8,
    keep_original=false)
```

For repeated applications to multiple validation series, fit once and reuse the fitted model:

```julia
model = fit_tvc(obs_series, raw_train_series; scales=[365, 183, 92, 46, 23, 12, 6, 3, 2])
corrected_series = apply_tvc(model, raw_val_series)
```

`tvc` returns a series on the validation time axis. The leading warm-up segment, of length `sum(scales) - length(scales)`, is filled with `NaN` because the multiscale rolling decomposition needs prior samples before the corrected series is well-defined.

## Practical Workflow

For scenario construction, the usual order is:

1. subset the calibration period
2. regrid observations and simulations to a common grid
3. choose a correction method
4. validate the corrected output before computing derived indicators

See [Building Climate Scenarios](scenarios.md) for the full workflow.

## Validation Checklist

After correction, compare the corrected data against the observational reference over the calibration period.

Useful checks include:

- changes in mean and standard deviation
- annual maxima and minima
- wet-day fractions or dry-spell lengths
- map patterns at representative dates
- time-series slices at representative grid points

See [Validation and Diagnostics](validation.md) for more detailed guidance.
