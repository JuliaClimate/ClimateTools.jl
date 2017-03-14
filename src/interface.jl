function Base.vcat(A::ClimGrid, B::ClimGrid)
  axisArray = vcat(A.data, B.data)
  ClimGrid(axisArray, model = A.model, experiment = A.experiment, run = A.run, filename = A.filename, dataunits = A.dataunits, latunits = A.latunits, lonunits = A.lonunits, var = A.var, typeof = A.typeof)
end
# TODO : Verify in Base.vcat(A::ClimGrid, B::ClimGrid) for time consistency (e.g. no same timestamp)
# TODO : Add methods for addition, subtraction, multiplication of ClimGrid types

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
    return C.var
  elseif i == 10
    return C.typeof
  else
    throw(error("You can't index like that!"))
  end
end

Base.linearindexing{T<:ClimGrid}(::Type{T}) = Base.LinearFast()
Base.size(::ClimGrid) = (10,)
Base.length(C::ClimGrid) = prod(size(C))
Base.endof(C::ClimGrid) = length(C)

function Base.show(C::ClimGrid)
  out = Dict()
  out["Data"] = string("Data array of size ", size(C[1].data))
  out["Model"] = C[3]
  out["Experiment"] = C[4]
  out["Run"] = C[5]
  out["Data Units"] = C[2]
  out["Latitude Units"] = C[7]
  out["Longitude Units"] = C[6]
  out["filename"] = C[8]
  out["var"] = C[9]
  out["typeof raw variable"] = C[10]

  return out
end
