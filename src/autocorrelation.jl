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
    indims = InDims("time")
    outdims = OutDims(Dim{:lags}(collect(1:lags)))    

    return mapCube(ClimateTools.autocorrelation, ds, indims=indims, outdims=outdims, lags=lags, nthreads=Threads.nthreads())


end

function hurst(xout, xin; k::Int=20)
    d = LongMemory.ClassicEstimators.rescaled_range_est(Array(xin), k=k)
    xout .= d + 0.5
    # xout .= LongMemory.autocorrelation(Array(xin), lags)   
end



"""
    hurst(ds; k)

Calculate the Hurst exponent for a time series of dataset `ds`. The Hurst exponent is a measure of the long-term memory of time series data. It can be used to classify time series as:

- **0.5 < H < 1**: The time series is persistent, meaning that high values are likely to be followed by high values, and low values by low values.
- **H = 0.5**: The time series is a random walk, indicating no correlation between values.
- **0 < H < 0.5**: The time series is anti-persistent, meaning that high values are likely to be followed by low values, and vice versa.

# Arguments
- `ds`: YAXArray variable.

# Returns
- `H`: The Hurst exponent of the time series.


"""
function hurst(ds; k::Int=20)

    # Dimensions
    indims = InDims("time")
    # outdims = OutDims(Dim{:hurst}(1))    
    outdims = OutDims()    

    return mapCube(ClimateTools.hurst, ds, indims=indims, outdims=outdims, k=k, nthreads=Threads.nthreads())



end
