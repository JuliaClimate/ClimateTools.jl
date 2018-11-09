# Bias correction

Quantile-quantile mapping (Themeßl et al. 2012, Grenier et al. 2015) is provided with ClimateTools.jl through the function [`qqmap`](@ref).

```julia
qqmap(obs::ClimGrid, ref::ClimGrid, fut::ClimGrid; method::String="Additive", detrend::Bool=true, window::Int=15, rankn::Int=50, thresnan::Float64=0.1, keep_original::Bool=false, interp = Linear(), extrap = Flat())
```

More information can be found in these references.

Themeßl, Matthias Jakob, Andreas Gobiet, and Georg Heinrich. 2012. “Empirical-Statistical Downscaling and Error Correction of Regional Climate Models and Its Impact on the Climate Change Signal.” Climatic Change 112 (2). Springer: 449–68.

Grenier, Patrick, Ramón de Elía, and Diane Chaumont. 2015. “Chances of Short-Term Cooling Estimated from a Selection of CMIP5-Based Climate Scenarios during 2006-2035 over Canada.” Journal of Climate, January 2015. American Meteorological Society. doi:10.1175/JCLI-D-14-00224.1.
