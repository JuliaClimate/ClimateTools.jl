# References

This page collects the main literature and methodological references used in the ClimateTools documentation.

## Bias Correction

- Themeßl, M. J., Gobiet, A., and Leuprecht, A. Empirical-statistical downscaling and error correction of regional climate models and its impact on the climate change signal. Theoretical and Applied Climatology, 2012.
- Piani, C., Haerter, J. O., and Coppola, E. Statistical bias correction for daily precipitation in regional climate models over Europe. Theoretical and Applied Climatology, 2010.
- Roy, P., Rondeau-Genesse, G., Jalbert, J., and Fournier, E. Climate scenarios of extreme precipitation using a combination of parametric and non-parametric bias correction methods in the province of Québec. Canadian Water Resources Journal, 2023. DOI: 10.1080/07011784.2023.2220682.

These references motivate the quantile-mapping workflows exposed through `qqmap` and related functions.

The Roy et al. (2023) paper is the methodological reference for `biascorrect_extremes`, including the QQM-GPD framing, external GEV parameter workflow, and transition between bulk and tail correction.

## Time Variability Correction

- Shao, Y., Bishop, C. H., Hobeichi, S., Nishant, N., Abramowitz, G., and Sherwood, S. Time Variability Correction of CMIP6 Climate Change Projections. Journal of Advances in Modeling Earth Systems, 2024. DOI: 10.1029/2023MS003640.

This is the methodological reference for `tvc`, `fit_tvc`, and `apply_tvc`.

## Climate Indicators

The grouped index interface in ClimateTools overlaps conceptually with ETCCDI and xclim-style climate-indicator workflows, especially for threshold counts, spell-duration indices, and percentile-based indices.

- Bourgault, P., Huard, D., Smith, T. J., Logan, T., Aoun, A., Lavoie, J., et al. xclim: xarray-based climate data analytics. Journal of Open Source Software, 2023. DOI: 10.21105/joss.05415.

This is the primary upstream reference for the xclim-style indicator framework that informs ClimateTools' grouped climate-index documentation and parity-oriented test coverage.

## Ensembles and Robustness

ClimateTools now includes a first xclim-inspired ensemble layer for summary statistics, percentile aggregation, deterministic subset selection, and robustness diagnostics.

- Bourgault, P., Huard, D., Smith, T. J., Logan, T., Aoun, A., Lavoie, J., et al. xclim: xarray-based climate data analytics. Journal of Open Source Software, 2023. DOI: 10.21105/joss.05415.
- Cannon, A. J. Selecting GCM Scenarios that Span the Range of Changes in a Multimodel Ensemble: Application to CMIP5 Climate Extremes Indices. Journal of Climate, 2015. DOI: 10.1175/JCLI-D-14-00636.1.
- Katsavounidis, I., Jay Kuo, C.-C., and Zhang, Z. A new initialization technique for generalized Lloyd iteration. IEEE Signal Processing Letters, 1994.
- Tebaldi, C., Arblaster, J. M., and Knutti, R. Mapping model agreement on future climate projections. Geophysical Research Letters, 2011. DOI: 10.1029/2011GL049863.
- Knutti, R., and Sedlacek, J. Robustness and uncertainties in the new CMIP5 climate model projections. Nature Climate Change, 2013. DOI: 10.1038/nclimate1716.
- Intergovernmental Panel on Climate Change (IPCC). Atlas. In Climate Change 2021: The Physical Science Basis, Cambridge University Press, 2023. DOI: 10.1017/9781009157896.021.

These references motivate the xclim-compatible ensemble APIs in ClimateTools: `ensemble_mean_std_max_min`, `ensemble_percentiles`, `make_criteria`, `kkz_reduce_ensemble`, `robustness_fractions`, `robustness_categories`, and `robustness_coefficient`.

## Regridding and Scenario Construction

The ClimateTools documentation uses workflow-oriented documentation patterns inspired by the broader climate-services ecosystem, including xclim’s separation of tutorials, workflow guides, and API reference.
