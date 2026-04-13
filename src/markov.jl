function _fit_msmodel(xin; k_regime=2, intercept="switching")
    return MarSwitching.MSModel(Float64.(collect(skipmissing(xin))), k_regime, intercept=intercept)
end

function MSModel(xout, xin; k_regime=2, intercept="switching")
    xout .= _fit_msmodel(xin; k_regime=k_regime, intercept=intercept)
end


"""
    MSModel(ds; k_switching=2, intercept="switching")

Compute the autocorrelation of a dataset `ds` along the time dimension.

# Arguments
- `ds`: The input dataset.
- `lags`: The number of lags to compute autocorrelation for. Default is 30.

# Returns
The autocorrelation of the input dataset `ds` along the time dimension.

"""
function MSModel(ds; k_regime=2, intercept = "switching")
    axis_names = Tuple(Symbol.(name.(ds.axes)))
    time_position = findfirst(==(:time), axis_names)
    time_position === nothing && error("MSModel requires a cube with a :time dimension.")

    retained_positions = Tuple(i for i in eachindex(ds.axes) if i != time_position)
    output_axes = tuple((ds.axes[i] for i in retained_positions)..., Dim{:MSM}(1:1))

    permutation = (retained_positions..., time_position)
    permuted_data = PermutedDimsArray(ds.data, permutation)
    retained_shape = size(permuted_data)[1:end-1]

    if isempty(retained_shape)
        fitted = _fit_msmodel(vec(permuted_data); k_regime=k_regime, intercept=intercept)
        return YAXArray(output_axes, reshape([fitted], 1))
    end

    fitted_models = Array{MarSwitching.MSM{Float64}}(undef, retained_shape..., 1)

    for index in CartesianIndices(Tuple(retained_shape))
        fitted_models[Tuple(index)..., 1] = _fit_msmodel(view(permuted_data, Tuple(index)..., :); k_regime=k_regime, intercept=intercept)
    end

    return YAXArray(output_axes, fitted_models)

end
