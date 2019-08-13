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
    mean(W::WeatherStation)

Compute the mean of `WeatherStation` W
"""
mean(W::WeatherStation) = Statistics.mean(filter(!isnan, W[1]))

"""
    minimum(A::ClimGrid)

Compute the minimum value of `ClimGrid` A
"""
Base.minimum(A::ClimGrid) = minimum(A[1])

"""
    minimum(W::WeatherStation)

Compute the minimum value of `WeatherStation` W
"""
Base.minimum(W::WeatherStation) = minimum(filter(!isnan, W[1]))

"""
    maximum(A::ClimGrid)

Compute the maximum value of `ClimGrid` A
"""
Base.maximum(A::ClimGrid) = maximum(A[1])

"""
    maximum(W::WeatherStation)

Compute the maximum value of `WeatherStation` W
"""
Base.maximum(W::WeatherStation) = maximum(filter(!isnan, W[1]))

"""
    std(A::ClimGrid)

Compute the standard deviation of `ClimGrid` A
"""
std(A::ClimGrid) = Statistics.std(A[1])

"""
    std(W::WeatherStation)

Compute the standard deviation of `WeatherStation` W
"""
std(W::WeatherStation) = Statistics.std(filter(!isnan, W[1]))

"""
    var(A::ClimGrid)

Compute the variance of `ClimGrid` A
"""
var(A::ClimGrid) = Statistics.var(A[1])

"""
    var(W::WeatherStation)

Compute the variance of `WeatherStation` W
"""
var(W::WeatherStation) = Statistics.var(filter(!isnan, W[1]))


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

function getindex(W::WeatherStation,i::Int)
    if i == 1
        return W.data
    elseif i == 2
        return W.dataunits
    elseif i == 3
        return W.lon
    elseif i == 4
        return W.lat
    elseif i == 5
        return W.alt
    elseif i == 6
        return W.lonunits
    elseif i == 7
        return W.latunits
    elseif i == 8
        return W.altunits
    elseif i == 9
        return W.stationID
    elseif i == 10
        return W.stationName
    elseif i == 11
        return W.filename
    elseif i == 12
        return W.variable
    elseif i == 13
        return W.typeofvar
    elseif i == 14
        return W.varattribs
    elseif i == 15
        return W.typeofcal
    elseif i == 16
        return W.timeattrib
    elseif i == 17
        return W.globalattribs
    else
        throw(error("You can't index like that!"))
    end
end

function getindex(W::WeatherNetwork,i::Int)
    return W.data[i]
end

function getindex(W::WeatherNetwork,s::String)
    idx = findall(x->x==s, W.stationID)[1]
    return W.data[idx]
end

# Base.IndexStyle{T<:ClimGrid}(::Type{T}) = Base.IndexLinear()
Base.length(C::ClimGrid) = length(fieldnames(typeof(C)))
Base.length(W::WeatherNetwork) = length(W.data)
Base.length(W::WeatherStation) = length(W.data)
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

Base.show(io::IO, ::MIME"text/plain", W::WeatherStation) = print(io, "WeatherStation struct with data:\n   ", summary(W[1]), "\n",
"Station ID: ", W[9], "\n",
"Station name: ", W[10], "\n",
"Variable: ", W[12], "\n",
"Data units: ", W[2], "\n",
"Global attributes: ", summary(W[17]), "\n",
"Filename: ", W[11])

Base.show(io::IO, ::MIME"text/plain", W::WeatherNetwork) = print(io, "WeatherNetwork struct with data:\n   ", summary(W[1]), "\n",
"Station IDs: ", W.stationID, "\n",
"Variable: ", W[1][12], "\n",
"Data units: ", W[1][2], "\n",
"Global attributes: ", summary(W[1][17]), "\n")

Base.show(io::IO, ::MIME"text/plain", ITP::TransferFunction) = print(io, "TransferFunction type with fields *itp*, *method* and *detrend*", "\n",
    "Interpolation array: ", size(ITP.itp), " transfer functions", "\n",
    "Method: ", ITP.method, "\n",
    "Detrended: ", ITP.detrend)
