# Bias correction

Quantile-quantile mapping (Themeßl et al. 2012, Grenier et al. 2015) is provided with ClimateTools.jl through the function [`qqmap`](@ref).

```julia
qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method::String="Additive", detrend::Bool=true, window::Int=15, rankn::Int=50, thresnan::Float64=0.1, keep_original::Bool=false, interp = Linear(), extrap = Flat())
```

Quantile-quantile mapping is a computational costly function. An alternative is to build a single transfer function over a domain.

```julia
qqmaptf(obs::ClimGrid, ref::ClimGrid; partition::Float64 = 1.0, method::String="Additive", window::Int64=15, rankn::Int64=50, interp = Linear(), extrap = Flat())
```

where `partition` is the fraction of grid points that will be considered in the construction of the transfer function. `qmaptf` returns `ITP` an Interpolation type. `ITP` then can be provided to the following `qqmap` method, which apply the transfer function `ITP` to values in `ClimGrid` `fut`.

```julia
qqmap(fut::ClimGrid, ITP; method::String="Additive")
```

More information can be found in these references.

Themeßl, Matthias Jakob, Andreas Gobiet, and Georg Heinrich. 2012. “Empirical-Statistical Downscaling and Error Correction of Regional Climate Models and Its Impact on the Climate Change Signal.” Climatic Change 112 (2). Springer: 449–68.

Grenier, Patrick, Ramón de Elía, and Diane Chaumont. 2015. “Chances of Short-Term Cooling Estimated from a Selection of CMIP5-Based Climate Scenarios during 2006-2035 over Canada.” Journal of Climate, January 2015. American Meteorological Society. doi:10.1175/JCLI-D-14-00224.1.
