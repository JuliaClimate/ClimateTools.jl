function estimate_prob(xout, xin1, xin2; da_cqt_quantile)
    
    # v_cqt = convert.(Float64,collect(da_cqt[31,31,:].data)) # the grid where we want to evaluate the interpolated value (ex. coords variable below)
    # v_eqt = collect(da_eqt[:,31,31].data) # the x-coordinate of the initial data
    
    # v_cqt_quantile = collect(da_cqt.quantile); # The y-coordinates of the data points, same length as `xp`.
    
    itp = linear_interpolation(Interpolations.deduplicate_knots!(xin2), da_cqt_quantile, extrapolation_bc=Flat());
    xout .= @. itp(xin1)
    
end


function estimate_prob(da_cqt::YAXArray, da_eq::YAXArray)

    Quantiles = collect(da_cqt.quantile)

    # Dimensions
    indims_da_cqt = InDims("quantile")
    indims_da_eqt = InDims("Quantiles")
    outdims = OutDims(Dim{:Probability}(Quantiles))

    return mapCube(estimate_prob, (da_cqt, da_eq), indims=(indims_da_cqt, indims_da_eqt), outdims=outdims, da_cqt_quantile=Quantiles)
    
end