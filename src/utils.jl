


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
    
    # On vérifie s'il y a un overlap sur l'année (e.g. on aimerait novembre, décember et janvier)
    if month2 < month1 # overlap annuel
        newcube1 = cube[time=Where(x -> month(x) >= month1)]
        newcube2 = cube[time=Where(x -> month(x) <= month2)]

        # newcube1 = renameaxis!(newcube1, :time=>:Ti)
        # newcube2 = renameaxis!(newcube2, :time=>:Ti)
        
        
        newlookupvals = vcat(collect(newcube1.time), collect(newcube2.time))
        # newlookupvals = vcat(collect(newcube1.Ti), collect(newcube2.Ti))

        # newcube = cat(newcube1, newcube2; dims=DimensionalData.Dimensions.Ti(newlookupvals))
        newcube = cat(newcube1, newcube2; dims=Dim{:time}(newlookupvals))
        # Dim{:time}(newlookupvals)


        # newcube = renameaxis!(newcube2, :Ti=>:time)
        
    else
        newcube = cube[time=Where(x -> month(x) <= month2 .&& month(x) >= month1)]
    end
    
    return newcube    
    
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
