function MSModel(xout, xin; k_regime=2, intercept="switching")    
    xout .= MarSwitching.MSModel(Float64.(xin), k_regime, intercept=intercept)
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
    return _xmap_call(
        Climat.MSModel,
        ds;
        reduced_dims=:time,
        output_axes=(Dim{:MSM}(1),),
        function_kwargs=(k_regime=k_regime, intercept=intercept),
    )

end
