


function dates_builder_year(x)
    out = Year[]
    for i in eachindex(x)
        push!(out, Year(x[i][1]))
    end

    return out
end

function dates_builder_yearmonth(x)
    out = Date[]
    for i in eachindex(x)
        push!(out, Date(x[i][1], x[i][2]))
    end

    return out
end

function dates_builder_yearmonth_hardcode(x, imois)
    out = Date[]
    for i in eachindex(x)
        push!(out, Date(x[i][1], imois))
    end

    return out
end



function dates_builder_yearmonthday(x)
    out = Date[]
    for i in eachindex(x)
        push!(out, Date(x[i][1], x[i][2], x[i][3]))
    end

    return out
end

function dates_builder_yearmonthday_hardcode(x; imois=1, iday=1)
    out = Date[]
    for i in eachindex(x)
        push!(out, Date(x[i][1], imois))
    end

    return out
end

function subsample(cube::YAXArray; month1=1, month2=12)
    
    # On vérifie s'il y a un overlap sur l'année (e.g. on aimerait novembre, décember et janvier)
    if month2 < month1 # overlap annuel
        @error "Not implemented. consider something along the lines of: newcube = cube[Ti=Where(x -> month(x) != month2)]"
        newdim = dim[Dim{:Ti}(Where(x -> month(x) <= month2 .| month(x) >= month1))]
        # index = (Dates.month.(datevecin) .<= endmonth) .&  (Dates.month.(datevecin) .>= startmonth)
    else
        newcube = cube[Ti=Where(x -> month(x) <= month2 .& month(x) >= month1)]        
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
