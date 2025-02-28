"""
    ERA5Land_dailysum(xout, xin; index_list = time_to_index)

Compute the daily sum of ERA5 land data. **used internally**

# Arguments
- `xout::AbstractArray`: Output array to store the computed daily sums.
- `xin::AbstractArray`: Input array containing ERA5 land data.
- `index_list::Function`: (optional) Function to convert time indices.

# Details
- The function computes the daily sum of ERA5 land data by doing first a diff operator and then summing the non-missing values in each day's data.
- Missing values are represented as `missing` in the input array and are replaced with `NaN` before computation.

"""
function ERA5Land_dailysum(xout, xin; index_list)
    
    xout .= NaN
    if !all(ismissing, xin)
        for i in eachindex(index_list)
            data = view(xin, index_list[i])
            # data = xin[index_list[i]]
            if !all(ismissing, data)
                diffxin = Base.diff(replace(data, missing => NaN))
                xout[i] = sum(diffxin[.!isnan.(diffxin)])
            end
        end        
    end    
end 

"""
    ERA5Land_dailysum(cube::YAXArray; keep_1stday=false, kwargs...)

Process ERA5-Land data by converting it to daily resolution. The function applies a diff function and then a sum for daily accumulated values.

# Arguments
- `cube::YAXArray`: The input data cube.
- `keep_1stday::Bool`: Whether to keep the first day of the data. Default is `false`.
- `kwargs...`: Additional keyword arguments.

# Returns
- `cube::YAXArray`: The processed data cube.

"""
function ERA5Land_dailysum(cube::YAXArray; keep_1stday=false, kwargs...)    
    
    # On veut des données quotidiennes pour ERA5-Land
    time_to_index = (cube.time .- Hour(1))
    time_index = yearmonthday.(time_to_index)
    new_dates = unique(time_index)
    
    if keep_1stday
        firstidx = 1        
    else
        firstidx = 2
    end
    index_in_cube = [findall(==(i), time_index) for i in unique(time_index)][firstidx:end]    
    
    # Dimensions
    indims = InDims("time")    
    outdims = OutDims(Dim{:time}(dates_builder_yearmonthday(new_dates[firstidx:end])))
    
    return mapCube(ERA5Land_dailysum, cube, indims=indims, outdims=outdims, index_list=index_in_cube)
    
end
