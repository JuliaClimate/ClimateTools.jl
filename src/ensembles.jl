"""
    ensemble_stats(xout, xin)

Compute the daily maximum, mean, and minimum values from the input array `xin` and store the results in the output array `xout`.

# Arguments
- `xout::AbstractArray`: The output array where the computed values will be stored.
- `xin::AbstractArray`: The input array from which the values will be computed.

# Details
- The function replaces missing values in `xin` with `NaN` before computing the maximum, mean, and minimum values.
- The computed values are stored in `xout` as an array, where each row corresponds to a different index in `index_list`.

"""
function ensemble_stats(xout, xin)
    xout .= NaN
    if !all(ismissing, xin)
        xin_missing = xin[.!ismissing.(xin)]
        xin_nan = xin_missing[.!isnan.(xin_missing)]
        xout[1] = Statistics.maximum(xin_nan)
        xout[2] = Statistics.mean(xin_nan)
        xout[3] = Statistics.minimum(xin_nan)
    end    
end

"""
        ensemble_stats(cube::YAXArray; dim="time")

Compute the maximum, mean, and minimum values from a YAXArray Cube along dimension `dim` (default to :time).

# Arguments
- `cube::YAXArray`: A Cube of data
- `dim::String`: The dimension along which to compute the statistics. Default is "Ti".

# Returns
- A Cube with maximum, mean and minimum values. The output cube has the following additional dimension:    
    - `stats`: max, mean and min   


"""
function ensemble_stats(cube::YAXArray; dim="time")    

      
    # Dimensions
    indims = InDims(dim)    
    # outdims = OutDims(RangeAxis("time", dates_builder_yearmonthday(new_dates)), CategoricalAxis("stats",["max","mean","min"]))
    outdims = OutDims(Dim{:stats}(["max","mean","min"]))
    
    mapCube(ensemble_stats, cube, indims=indims, outdims=outdims)
end
