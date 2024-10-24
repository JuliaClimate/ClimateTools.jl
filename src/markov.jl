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

    # Dimensions
    indims = InDims("Ti")
    outdims = OutDims(Dim{:MSM}(1))
    # outdims = OutDims()

    return mapCube(Climat.MSModel, ds, indims=indims, outdims=outdims, k_regime=k_regime, intercept=intercept, nthreads=Threads.nthreads())

end
