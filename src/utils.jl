


function dates_builder_year(x)
    out = Year[]
    for i in eachindex(x)
        push!(out, Year(x[i][1]))
    end

    return out
end

function dates_builder_yearmonth(x)
    out = DateTime[]
    for i in eachindex(x)
        push!(out, DateTime(x[i][1], x[i][2]))
    end

    return out
end

function dates_builder_yearmonth_hardcode(x, imois)
    out = DateTime[]
    for i in eachindex(x)
        push!(out, DateTime(x[i][1], imois))
    end

    return out
end

function dates_builder_monthly_resample(x)
    out = DateTime[]
    for i in eachindex(x)
        push!(out, DateTime(x[i][1], x[i][2], 15))
    end

    return out
end

function dates_builder_yearmonthday(x)
    out = DateTime[]
    for i in eachindex(x)
        push!(out, DateTime(x[i][1], x[i][2], x[i][3]))
    end

    return out
end

function dates_builder_yearmonthday_hardcode(x; imois=1, iday=1)
    out = DateTime[]
    for i in eachindex(x)
        push!(out, DateTime(x[i][1], imois))
    end

    return out
end

function subsample(cube::YAXArray; month1=1, month2=12)
    predicate = if month2 < month1
        x -> month(x) >= month1 || month(x) <= month2
    else
        x -> month1 <= month(x) <= month2
    end

    return cube[time=Where(predicate)]
end

"""
    m2mm(cube::YAXArray, kwargs...)

Converts the values in the input `cube` from meters to millimeters.

# Arguments
- `cube::YAXArray`: The input cube containing values in meters.
- `kwargs...`: Additional keyword arguments.

# Returns
- `cubeout`: The output cube with values converted to millimeters.

"""
function m2mm(cube::YAXArray, kwargs...)
    
    cubeout = map(cube) do x
        x * 1000.0
    end

    return cubeout
    
end
