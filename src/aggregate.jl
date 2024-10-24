# function Base.sum(xout, xin)    
    
#     xout .=  Base.sum(xin)
   
# end

# function Base.sum(cube::YAXArray, kwargs...)    

#     indims = InDims("time")    
#     outdims = OutDims()
    
#     mapCube(sum, cube, indims=indims, outdims=outdims)
# end

"""
    daily_max_mean_min(xout, xin)

Compute the daily maximum, mean, and minimum values from the input array `xin` and store the results in the output array `xout`.

# Arguments
- `xout::AbstractArray`: The output array where the computed values will be stored.
- `xin::AbstractArray`: The input array from which the values will be computed.

# Details
- The function replaces missing values in `xin` with `NaN` before computing the maximum, mean, and minimum values.
- The computed values are stored in `xout` as an array, where each row corresponds to a different index in `index_list`.

"""
function daily_max_mean_min(xout, xin; index_list)
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            data = view(xin, index_list[i])
            if !all(ismissing, data)
                xout[i,1] = maximum(skipmissing(data))
                xout[i,2] = mean(skipmissing(data))
                xout[i,3] = minimum(skipmissing(data))
            end
        end        
    end    
end

"""
        daily_max_mean_min(data::YAXArray)

Compute the daily maximum, mean, and minimum values from a YAXArray Cube.

# Arguments
- `data::YAXArray`: A Cube of data

# Returns
- A Cube with daily maximum, mean and minimum values. The output cube has the following dimensions:
    - `time`: Time dimension
    - `stats`: max, mean and min   


"""
function daily_max_mean_min(cube::YAXArray, kwargs...)    

    time_to_index = cube.time
    time_index = yearmonthday.(time_to_index)
    new_dates = unique(time_index)
    index_in_cube = [findall(==(i), time_index) for i in unique(time_index)]   
    
    # Dimensions
    indims = InDims("Ti")    
    # outdims = OutDims(RangeAxis("time", dates_builder_yearmonthday(new_dates)), CategoricalAxis("stats",["max","mean","min"]))
    outdims = OutDims(Dim{:Ti}(dates_builder_yearmonthday(new_dates)), CategoricalAxis("stats",["max","mean","min"]))
    
    mapCube(daily_max_mean_min, cube, indims=indims, outdims=outdims, index_list=index_in_cube)
end


"""
    percentiles(cube::YAXArray, quantiles=[0.1, 0.5, 0.9])

Compute the percentiles of a YAXArray cube along the "number" dimension (represent member for a given seasonal prediction simulation) and return a new cube with the computed percentiles along the "time", "longitude", and "latitude" dimensions.

# Arguments
- `cube::YAXArray`: The input YAXArray cube.
- `quantiles::Vector{Float64}`: The quantiles to compute. Default is `[0.1, 0.5, 0.9]`.

# Returns
A new YAXArray cube with the computed percentiles along the "time", "longitude", and "latitude" dimensions.

"""
function percentiles(cube::YAXArray, quantiles=[0.1, 0.5, 0.9];lonname="longitude", latname="latitude")    
    
    indims = InDims("number")    
    outdims = OutDims("Ti", lonname, latname)
    
    mapCube(percentiles, cube, indims=indims, outdims=outdims, quantiles=quantiles)
end

"""
    diff(xout, xin)

Compute the difference between consecutive elements of `xin` and store the result in `xout`. A 0.0 is added at the beginning of the output array.

# Arguments
- `xout`: An array to store the difference values.
- `xin`: An array of values.

"""
function diff(xout, xin)   
        
    xout .=  vcat(0.0, Base.diff(xin[:]))    
   
end

"""
    diff(cube::YAXArray, kwargs...)

Compute the difference between consecutive time steps in the given `cube`.

# Arguments
- `cube::YAXArray`: The input data cube.
- `kwargs...`: Additional keyword arguments.

# Returns
- A new `YAXArray` object containing the difference between consecutive time steps.

"""
function diff(cube::YAXArray)    


    indims = InDims("Ti")    
    outdims = OutDims("Ti")
    
    mapCube(diff, cube, indims=indims, outdims=outdims)
end

"""
    cumsum(xout, xin)

Compute the cumulative sum of the elements in `xin` and store the result in `xout`.

# Arguments
- `xout`: An array to store the cumulative sum.
- `xin`: An array containing the elements to be summed.

"""
function cumsum(xout, xin)    
    
    xout .=  cumsum(xin)
   
end

"""
    cumsum(cube::YAXArray, kwargs...)

Compute the cumulative sum along the "time" dimension of the input `cube`.

# Arguments
- `cube`: The input YAXArray cube.
- `kwargs`: Additional keyword arguments.

# Returns
A new YAXArray cube with the cumulative sum computed along the "time" dimension.

"""
function cumsum(cube::YAXArray, kwargs...)    

    indims = InDims("Ti")    
    outdims = OutDims("Ti")
    
    mapCube(cumsum, cube, indims=indims, outdims=outdims)
end



"""
    daily_fct(xout, xin; fct::Function, index_list = time_to_index)

Apply a function `fct` to aggregate data from `xin` into `xout` on a daily basis **used internally**

# Arguments
- `xout`: Output array to store the aggregated data.
- `xin`: Input array containing the data to be aggregated.
- `fct`: Function to be applied for aggregation.
- `index_list`: (optional) Function to convert time to index.

"""
function daily_fct(xout, xin; fct::Function, index_list = time_to_index)
    
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            data = view(xin, index_list[i])
            if !all(ismissing, data)
                xout[i] = fct(skipmissing(data))
            end
        end
        
    end
    
end 

"""
    daily_fct(cube::YAXArray; fct::Function=mean, shifthour=0, kwargs...)

Apply a function `fct` to process the data in `cube` on a daily basis. Assumes that the Cube have sub-daily values.

# Arguments
- `cube::YAXArray`: The input data cube.
- `fct::Function`: The function to apply for processus sub-daily values. Default is `mean`.
- `shifthour::Int`: The number of hours to shift the time values. Default is 0.
- `kwargs...`: Additional keyword arguments to be passed to the aggregation function.

# Returns
A new `YAXArray` with aggregated data on a daily basis.

"""
function daily_fct(cube::YAXArray; fct::Function=mean, shifthour=0, kwargs...)
    time_to_index = cube.Ti .+ Dates.Hour(shifthour)
    time_index = yearmonthday.(time_to_index)
    new_dates = unique(time_index)
    index_in_cube = [findall(==(i), time_index) for i in unique(time_index)]   
    
    # Dimensions
    indims = InDims("Ti")        
    outdims = OutDims(Dim{:Ti}(dates_builder_yearmonthday(new_dates)))
    
    mapCube(daily_fct, cube, indims=indims, outdims=outdims, fct=fct, index_list=index_in_cube)
end


function yearly_clim(xout, xin; fct::Function=mean, index_list = time_to_index)
    
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            data = view(xin, index_list[i])
            if !all(ismissing, data)
                xout[i] = fct(skipmissing(data))
            end
        end
        
    end
    
end 
"""
    yearly_clim(cube::YAXArray; fct::Function=mean, kwargs...)

Compute the yearly climatology of a YAXArray.

# Arguments
- `cube::YAXArray`: The input YAXArray containing the data.
- `fct::Function`: The function to be applied to each yearly subset of the data. Default is `mean`.
- `kwargs...`: Additional keyword arguments to be passed to the function. **not yet implemented**

# Returns
A new YAXArray containing the yearly climatology.
"""

function yearly_clim(cube::YAXArray; fct::Function=mean, kwargs...)
    time_to_index = cube.Ti;
    time_index = year.(time_to_index);
    new_dates = unique(time_index);
    index_in_cube = [findall(==(i), time_index) for i in unique(time_index)]

    # Dimensions
    indims = InDims("Ti")        
    outdims = OutDims(Dim{:Ti}(dates_builder_yearmonthday_hardcode(new_dates, imois=7, iday=1)))

    mapCube(yearly_clim, cube, fct=fct, indims=indims, outdims=outdims, index_list=index_in_cube)
end

"""
    function daily_max(xout, xin; index_list = time_to_index)

xout = expected output
xin = input that I get from the Cube
index_list = 
"""
function daily_max(xout, xin; index_list = time_to_index)
    
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            data = view(xin, index_list[i])
            if !all(ismissing, data)
                xout[i] = maximum(skipmissing(data))
            end
        end
        
    end
    
end 

function daily_max(cube::YAXArray, kwargs...)
    time_to_index = cube.time;
    time_index = yearmonthday.(time_to_index);
    new_dates = unique(time_index);
    index_in_cube = [findall(==(i), time_index) for i in unique(time_index)]   
    
    # Dimensions
    indims = InDims("Ti")        
    outdims = OutDims(Dim{:Ti}(dates_builder_yearmonthday(new_dates)))
    
    mapCube(daily_max, cube, indims=indims, outdims=outdims, index_list=index_in_cube)
end
