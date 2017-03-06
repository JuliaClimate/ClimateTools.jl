function Base.vcat(A::ClimGrid, B::ClimGrid)
  axisArray = vcat(A.data, B.data)
  ClimGrid(axisArray, A.model, A.experiment, A.run, A.filename, A.dataunits, A.latunits, A.lonunits)
end
