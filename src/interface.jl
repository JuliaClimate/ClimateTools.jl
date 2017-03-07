function Base.vcat(A::ClimGrid, B::ClimGrid)
  axisArray = vcat(A.data, B.data)
  ClimGrid(axisArray, A.model, A.experiment, A.run, A.filename, A.dataunits, A.latunits, A.lonunits)
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
    return C.var
  else
    throw(error("You can't index like that!"))
  end
end

Base.linearindexing{T<:ClimGrid}(::Type{T}) = Base.LinearFast()
Base.size(C::ClimGrid) = (8,)
Base.length(C::ClimGrid) = prod(size(C))
Base.endof(C::ClimGrid) = length(C)

function Base.show(C::ClimGrid)
  out = Dict()
  # out["Data"] = C[1]
  out["Model"] = C[3]
  out["Experiment"] = C[4]
  out["Run"] = C[5]
  out["Data Units"] = C[2]
  out["Latitude Units"] = C[7]
  out["Longitude Units"] = C[6]
  out["Filename"] = C[8]

  return out

  # s1 = string("      Data = AxisArray with size ", size(C[1]), " and axis names ", axisnames(C[1]))
  # s2 = string("     model = ", C[3])
  # s3 = string("experiment = ", C[4])
  # s4 = string("       run = ", C[5])
  # s5 = string(" dataunits = ", C[2])
  # s6 = string("Lat. units = ", C[7])
  # s7 = string("Lon. units = ", C[6])
  # s8 = string("  Filename = ", C[8])
  # return [println(s1), println(s2), println(s3), println(s4), println(s5), println(s6), println(s7), println(s8)];

  # return vcat(s1, s2, s3, s4, s5, s6, s7, s8)
end
