function autocorrelation(xout, xin; lags=30)    
    xout .= LongMemory.autocorrelation(Array(xin), lags)   
end


"""
    autocorrelation(ds; lags=30)

Compute the autocorrelation of a dataset `ds` along the time dimension.

# Arguments
- `ds`: The input dataset.
- `lags`: The number of lags to compute autocorrelation for. Default is 30.

# Returns
The autocorrelation of the input dataset `ds` along the time dimension.

"""
function autocorrelation(ds; lags=30)

    # Dimensions
    indims = InDims("Ti")
    outdims = OutDims(Dim{:lags}(collect(1:lags)))    

    return mapCube(ClimateTools.autocorrelation, ds, indims=indims, outdims=outdims, lags=lags, nthreads=Threads.nthreads())

end
