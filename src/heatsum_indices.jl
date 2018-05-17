"""
    daysabove10(C::ClimGrid)

Annual number of days with temperature >= 10 Celsius. This function returns a ClimGrid.
"""

# TODO add check for units. important for threshold indices

function daysabove10(C::ClimGrid)
  @argcheck C[9] == "tas"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (size(C.data, 1), size(C.data, 2), length(numYears)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 10, +, view(dataout, :, :, i:i), view(datain, :, :, idx))
  end

  # Build output AxisArray
  FD = buildarray(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits="Days", latunits=C.latunits, lonunits=C.lonunits, variable="daysabove10", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

# """
#     daysabove10(data::AbstractArray, timevector::StepRange{Date,Base.Dates.Day})
#
# Annual number with temperature >= 10 Celsius. This function returns a boolean vector. *true* if the data is higher or equal to 1 and *false* otherwise.
# """
#
# function daysabove10(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
#   years    = Dates.year.(timeV)
#   numYears = unique(years)
#   FD       = zeros(Int64, (length(numYears)))
#
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t >= 10, +, view(FD, i:i), view(data, idx))
#   end
#   return FD
# end
#
# function daysabove10(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
#   years    = Dates.year.(timeV)
#   numYears = unique(years)
#   FD       = zeros(Int64, (size(data, 1), length(numYears)))
#
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t >= 10, +, view(FD, :, i:i), view(data, :, :, idx))
#   end
#   return FD
# end
#
# function daysabove10(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
#   years    = Dates.year.(timeV)
#   numYears = unique(years)
#   FD       = zeros(Int64, (size(data, 2), size(data, 3), length(numYears)))
#
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t >= 10, +, view(FD, :, :, i:i), view(data, :, :, idx))
#   end
#   return FD
# end
