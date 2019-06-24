function buildarrayinterface(axisArraytmp, A)
    latsymbol = Symbol(A.dimension_dict["lat"])
    lonsymbol = Symbol(A.dimension_dict["lon"])
    if ndims(axisArraytmp) == 2
        axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}][:]), Axis{latsymbol}(A[1][Axis{latsymbol}][:]))
    elseif ndims(axisArraytmp) == 3
        axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}][:]), Axis{latsymbol}(A[1][Axis{latsymbol}][:]), Axis{:time}(A[1][Axis{:time}][:]))
    elseif ndims(axisArraytmp) == 4
        axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}][:]), Axis{latsymbol}(A[1][Axis{latsymbol}][:]), Axis{:plev}(A[1][Axis{:plev}][:]), Axis{:time}(A[1][Axis{:time}][:]))
    end
    return axisArray
end

function getsymbols(C::ClimGrid)
    latsymbol = Symbol(C.dimension_dict["lat"])
    lonsymbol = Symbol(C.dimension_dict["lon"])

    return latsymbol, lonsymbol
end

"""
    merge(A::ClimGrid, B::ClimGrid)

Combines two ClimGrid. Based on the AxisArrays method.
"""
function Base.merge(A::ClimGrid, B::ClimGrid)
    axisArray = merge(A.data, B.data)
    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=A.experiment, run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:+(A::ClimGrid, B::ClimGrid)
    axisArraytmp = A.data .+ B.data

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " + ", B.experiment), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:+(A::ClimGrid, k)
    axisArraytmp = A.data .+ k

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " + ", k), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:-(A::ClimGrid, B::ClimGrid)

    axisArraytmp = A.data .- B.data

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " - ", B.experiment), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:-(A::ClimGrid, k)
    axisArraytmp = A.data .- k

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " - ", k), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:*(A::ClimGrid, B::ClimGrid)

    axisArraytmp = A.data .* B.data

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " * ", B.experiment), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:*(A::ClimGrid, k)
    axisArraytmp = A.data .* k

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " * ", k), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:/(A::ClimGrid, B::ClimGrid)
    axisArraytmp = A.data ./ B.data

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " / ", B.experiment), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:/(A::ClimGrid, k)
    axisArraytmp = A.data ./ k

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " / ", k), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

"""
    mean(A::ClimGrid)

Compute the mean of `ClimGrid` A
"""
mean(A::ClimGrid) = Statistics.mean(A[1])

"""
    minimum(A::ClimGrid)

Compute the minimum value of `ClimGrid` A
"""
Base.minimum(A::ClimGrid) = minimum(A[1])

"""
    maximum(A::ClimGrid)

Compute the maximum value of `ClimGrid` A
"""
Base.maximum(A::ClimGrid) = maximum(A[1])

"""
    std(A::ClimGrid)

Compute the standard deviation of `ClimGrid` A
"""
std(A::ClimGrid) = Statistics.std(A[1])

"""
    var(A::ClimGrid)

Compute the variance of `ClimGrid` A
"""
var(A::ClimGrid) = Statistics.var(A[1])


function getindex(C::ClimGrid,i::Int)
    if i == 1
        return C.data
    elseif i == 2
        return C.dataunits
    elseif i == 3
        return C.model
    elseif i == 4
        return C.experiment
    elseif i == 5
        return C.run
    elseif i == 6
        return C.lonunits
    elseif i == 7
        return C.latunits
    elseif i == 8
        return C.filename
    elseif i == 9
        return C.variable
    elseif i == 10
        return C.typeofvar
    elseif i == 11
        return C.typeofcal
    elseif i == 12
        return C.globalattribs
    else
        throw(error("You can't index like that!"))
    end
end

function getindex(C::WeatherStation,i::Int)
    if i == 1
        return C.data
    elseif i == 2
        return C.dataunits
    elseif i == 3
        return C.lon
    elseif i == 4
        return C.lat
    elseif i == 5
        return C.alt
    elseif i == 6
        return C.lonunits
    elseif i == 7
        return C.latunits
    elseif i == 8
        return C.altunits
    elseif i == 9
        return C.stationID
    elseif i == 10
        return C.stationName
    elseif i == 11
        return C.filename
    elseif i == 12
        return C.variable
    elseif i == 13
        return C.typeofvar
    elseif i == 14
        return C.varattribs
    elseif i == 15
        return C.typeofcal
    elseif i == 16
        return C.timeattrib
    elseif i == 17
        return C.globalattribs
    else
        throw(error("You can't index like that!"))
    end
end

function getindex(C::WeatherNetwork,i::Int)
    return C.data[i]
end

function getindex(C::WeatherNetwork,s::String)
    idx = findall(x->x==s, C.stationID)[1]
    return C.data[idx]
end

# Base.IndexStyle{T<:ClimGrid}(::Type{T}) = Base.IndexLinear()
Base.length(C::ClimGrid) = length(fieldnames(typeof(C)))
Base.size(C::ClimGrid) = (length(C),)
Base.size(C::ClimGrid, n::Int) = n==1 ? length(C) : error("Only dimension 1 has a well-defined size.")
#Base.endof(C::ClimGrid) = length(C)
Base.ndims(::ClimGrid) = 1


Base.show(io::IO, ::MIME"text/plain", C::ClimGrid) = print(io, "ClimGrid struct with data:\n   ", summary(C[1]), "\n",
"Project: ", C.project, "\n",
"Institute: ", C.institute, "\n",
"Model: ", C[3], "\n",
"Experiment: ", C[4], "\n",
"Run: ", C[5], "\n",
"Variable: ", C[9], "\n",
# "Variable CF standard name: ", C.varattribs["standard_name"], "\n",
"Data units: ", C[2], "\n",
"Frequency: ", C.frequency, "\n",
"Global attributes: ", summary(C[12]), "\n",
"Filename: ", C[8])

Base.show(io::IO, ::MIME"text/plain", C::WeatherStation) = print(io, "WeatherStation struct with data:\n   ", summary(C[1]), "\n",
"Station ID: ", C[9], "\n",
"Station name: ", C[10], "\n",
"Variable: ", C[12], "\n",
"Data units: ", C[2], "\n",
"Global attributes: ", summary(C[17]), "\n",
"Filename: ", C[11])

Base.show(io::IO, ::MIME"text/plain", ITP::TransferFunction) = print(io, "TransferFunction type with fields *itp*, *method* and *detrend*", "\n",
    "Interpolation array: ", size(ITP.itp), " transfer functions", "\n",
    "Method: ", ITP.method, "\n",
    "Detrended: ", ITP.detrend)
