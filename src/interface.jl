"""
    vcat(A::ClimGrid, B::ClimGrid)

Combines two ClimGrid. Based on the AxisArrays method. Better way to do it would be to use the merge method.
"""

function Base.vcat(A::ClimGrid, B::ClimGrid)
    warn("Use merge function instead of vcat for dimensions consistency")
    axisArraytmp = vcat(A.data, B.data)
    # TODO add axis information in the creation of the AxisArray
    axisArray = AxisArray(axisArraytmp)
    ClimGrid(axisArray, model = A.model, experiment = A.experiment, run = A.run, filename = A.filename, dataunits = A.dataunits, latunits = A.latunits, lonunits = A.lonunits, variable = A.variable, typeofvar = A.typeofvar, typeofcal = A.typeofcal)
end

"""
    merge(A::ClimGrid, B::ClimGrid)

Combines two ClimGrid. Based on the AxisArrays method.
"""

function Base.merge(A::ClimGrid, B::ClimGrid)
    axisArray = merge(A.data, B.data)
    ClimGrid(axisArray, model = A.model, experiment = A.experiment, run = A.run, filename = A.filename, dataunits = A.dataunits, latunits = A.latunits, lonunits = A.lonunits, variable = A.variable, typeofvar = A.typeofvar, typeofcal = A.typeofcal)
end

function Base.:+(A::ClimGrid, B::ClimGrid)
    axisArraytmp = A.data + B.data

    axisArray = AxisArray(axisArraytmp, Axis{:time}(A[1][Axis{:time}][:]), Axis{:lon}(A[1][Axis{:lon}][:]), Axis{:lat}(A[1][Axis{:lat}][:]))

    ClimGrid(axisArray, model = A.model, experiment = A.experiment, run = A.run, filename = A.filename, dataunits = A.dataunits, latunits = A.latunits, lonunits = A.lonunits, variable = A.variable, typeofvar = A.typeofvar, typeofcal = A.typeofcal)
end

function Base.:-(A::ClimGrid, B::ClimGrid)
    axisArraytmp = A.data - B.data

    axisArray = AxisArray(axisArraytmp, Axis{:time}(A[1][Axis{:time}][:]), Axis{:lon}(A[1][Axis{:lon}][:]), Axis{:lat}(A[1][Axis{:lat}][:]))

    ClimGrid(axisArray, model = A.model, experiment = string(A.experiment, " minus ", B.experiment), run = A.run, filename = A.filename, dataunits = A.dataunits, latunits = A.latunits, lonunits = A.lonunits, variable = A.variable, typeofvar = A.typeofvar, typeofcal = A.typeofcal)
end

function Base.:*(A::ClimGrid, k)
    axisArraytmp = A.data * k

    axisArray = AxisArray(axisArraytmp, Axis{:time}(A[1][Axis{:time}][:]), Axis{:lon}(A[1][Axis{:lon}][:]), Axis{:lat}(A[1][Axis{:lat}][:]))

    ClimGrid(axisArray, model = A.model, experiment = string(A.experiment, " multiplied by ", k), run = A.run, filename = A.filename, dataunits = A.dataunits, latunits = A.latunits, lonunits = A.lonunits, variable = A.variable, typeofvar = A.typeofvar, typeofcal = A.typeofcal)
end

# function Base.:*(A::ClimGrid, B::ClimGrid)
#   axisArraytmp = A.data * B.data
#
#   axisArray = AxisArray(axisArraytmp, Axis{:time}(A[1][Axis{:time}][:]), Axis{:lon}(A[1][Axis{:lon}][:]), Axis{:lat}(A[1][Axis{:lat}][:]))
#
#   ClimGrid(axisArray, model = A.model, experiment = string(A.experiment, " multiplied by another ClimGrid"), run = A.run, filename = A.filename, dataunits = A.dataunits, latunits = A.latunits, lonunits = A.lonunits, variable = A.variable, typeofvar = A.typeofvar, typeofcal = A.typeofcal)
# end

function Base.:/(A::ClimGrid, k)
    axisArraytmp = A.data / k

    axisArray = AxisArray(axisArraytmp, Axis{:time}(A[1][Axis{:time}][:]), Axis{:lon}(A[1][Axis{:lon}][:]), Axis{:lat}(A[1][Axis{:lat}][:]))

    ClimGrid(axisArray, model = A.model, experiment = string(A.experiment, " divided by ", k), run = A.run, filename = A.filename, dataunits = A.dataunits, latunits = A.latunits, lonunits = A.lonunits, variable = A.variable, typeofvar = A.typeofvar, typeofcal = A.typeofcal)
end

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
    else
        throw(error("You can't index like that!"))
    end
end

Base.IndexStyle{T<:ClimGrid}(::Type{T}) = Base.IndexLinear()
Base.length(C::ClimGrid) = length(fieldnames(C))
Base.size(C::ClimGrid) = (length(C),)
Base.size(C::ClimGrid,n::Int) = n==1 ? length(C) : error("Only dimension 1 has a well-defined size.")
Base.endof(C::ClimGrid) = length(C)
Base.ndims(::ClimGrid) = 1



Base.show(io::IO, ::MIME"text/plain", C::ClimGrid) =
           print(io, "ClimGrid struct with data:\n   ", summary(C[1]), "\n", "Model: ", C[3], "\n", "Experiment: ", C[4], "\n", "Run: ", C[5], "\n", "Variable: ", C[9], "\n", "Data units: ", C[2], "\n", "Filename: ", C[8])
