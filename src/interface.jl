function buildarrayinterface(axisArraytmp, A)
    latsymbol = Symbol(A.dimension_dict["lat"])
    lonsymbol = Symbol(A.dimension_dict["lon"])
    axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}][:]), Axis{latsymbol}(A[1][Axis{latsymbol}][:]), Axis{:time}(A[1][Axis{:time}][:]))
    return axisArray
end

# """
#     vcat(A::ClimGrid, B::ClimGrid)
#
# Combines two ClimGrid. Based on the AxisArrays method. Better way to do it would be to use the merge method.
# """
#
# function Base.vcat(A::ClimGrid, B::ClimGrid)
#     warn("Use merge function instead of vcat for dimensions consistency")
#     axisArraytmp = vcat(A.data, B.data, 3)
#     # TODO add axis information in the creation of the AxisArray
#     axisArray = AxisArray(axisArraytmp)
#     ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, model=A.model, frequency=A.frequency, experiment=A.experiment, run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
# end

"""
    merge(A::ClimGrid, B::ClimGrid)

Combines two ClimGrid. Based on the AxisArrays method.
"""
function Base.merge(A::ClimGrid, B::ClimGrid)
    axisArray = merge(A.data, B.data)
    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, model=A.model, frequency=A.frequency, experiment=A.experiment, run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:+(A::ClimGrid, B::ClimGrid)
    axisArraytmp = A.data + B.data

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " + ", B.experiment), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:+(A::ClimGrid, k)
    axisArraytmp = A.data + k

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " + ", k), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:-(A::ClimGrid, B::ClimGrid)
    axisArraytmp = A.data - B.data

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " - ", B.experiment), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:-(A::ClimGrid, k)
    axisArraytmp = A.data - k

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " - ", k), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:*(A::ClimGrid, B::ClimGrid)
    axisArraytmp = A.data .* B.data

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " * ", B.experiment), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:*(A::ClimGrid, k)
    axisArraytmp = A.data .* k

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " * ", k), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:/(A::ClimGrid, B::ClimGrid)
    axisArraytmp = A.data ./ B.data

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " / ", B.experiment), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

function Base.:/(A::ClimGrid, k)
    axisArraytmp = A.data / k

    axisArray = buildarrayinterface(axisArraytmp, A)

    ClimGrid(axisArray, longrid=A.longrid, latgrid=A.latgrid, msk=A.msk, grid_mapping=A.grid_mapping, dimension_dict=A.dimension_dict, model=A.model, frequency=A.frequency, experiment=string(A.experiment, " / ", k), run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=A.latunits, lonunits=A.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)
end

"""
    mean(A::ClimGrid)

Compute the mean of `ClimGrid` A
"""
mean(A::ClimGrid) = mean(A[1])

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
std(A::ClimGrid) = std(A[1])

"""
    var(A::ClimGrid)

Compute the variance of `ClimGrid` A
"""
var(A::ClimGrid) = var(A[1])


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

# Base.IndexStyle{T<:ClimGrid}(::Type{T}) = Base.IndexLinear()
Base.length(C::ClimGrid) = length(fieldnames(C))
Base.size(C::ClimGrid) = (length(C),)
Base.size(C::ClimGrid, n::Int) = n==1 ? length(C) : error("Only dimension 1 has a well-defined size.")
Base.endof(C::ClimGrid) = length(C)
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

Base.show(io::IO, ::MIME"text/plain", ITP::TransferFunction) = print(io, "TransferFunction type with fields *itp*, *method* and *detrend*", "\n",
    "Interpolation array: ", size(ITP.itp), " transfer functions", "\n",
    "Method: ", ITP.method, "\n",
    "Detrended: ", ITP.detrend)


"""
    timeindex(timeVec, start_date, end_date, freq)

Return the index of time vector specified by start_date and end_date. Provide timestep "freq" to account for monthly timestep.
"""
function timeindex(timeV, datebeg::Tuple, dateend::Tuple, frequency)

    # Start Date
    if !isinf(datebeg[1])

        # Build DateTime type
        start_date = buildtimetype(datebeg)

        # Check
        @argcheck start_date >= timeV[1]

        if frequency != "mon"
            idxtimebeg = findfirst(timeV .== start_date)[1]
        elseif frequency == "mon"
            idxtimebeg = findfirst(timeV .== start_date)[1]
        end
    else
        idxtimebeg = 1

    end

    # End date
    if !isinf(dateend[1])

        # Build DateTime type
        end_date = buildtimetype(dateend)

        @argcheck end_date <= timeV[end]
        if frequency != "mon"
            idxtimeend = findlast(timeV .== end_date)[1]
        elseif frequency == "mon"
            idxtimeend = findlast(timeV .== Date(Dates.year(end_date), Dates.month(end_date), Dates.day(1)))[1]
        end
    else
        idxtimeend = length(timeV)
    end

    if !isinf(datebeg[1]) && !isinf(dateend[1])
        @argcheck start_date <= end_date
    end

    return idxtimebeg, idxtimeend
end

"""
    buildtimetype(datetuple)

Returns the adequate DateTime for temporal subsetting.
"""
function buildtimetype(datetuple)

    if length(datetuple) == 1
        dateout = DateTime(datetuple[1], 01, 01)
    elseif length(datetuple) == 2
        dateout = DateTime(datetuple[1], datetuple[2], 01)
    elseif length(datetuple) == 3
        dateout = DateTime(datetuple[1], datetuple[2], datetuple[3])
    elseif length(datetuple) == 4
        dateout = DateTime(datetuple[1], datetuple[2], datetuple[3], datetuple[4], 00, 00)
    elseif length(datetuple) == 5
        dateout = DateTime(datetuple[1], datetuple[2], datetuple[3], datetuple[4], datetuple[5], 00)
    elseif length(datetuple) == 6
        dateout = DateTime(datetuple[1], datetuple[2], datetuple[3], datetuple[4], datetuple[5], datetuple[6])
    end

    return dateout


end
