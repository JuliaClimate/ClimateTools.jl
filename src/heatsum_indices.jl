function buildarray(C::ClimateTools.ClimGrid, dataout, numYears)
    lonsymbol = Symbol(C.dimension_dict["lon"])
    latsymbol = Symbol(C.dimension_dict["lat"])
    FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{lonsymbol}(C[1][Axis{lonsymbol}][:]), Axis{latsymbol}(C[1][Axis{latsymbol}][:]))
    return FD
end

"""
    daysabove10(C::ClimGrid)

Annual number of days with temperature >= 10 Celsius. This function returns a ClimGrid.
"""

# TODO add check for units. important for threshold indices

function daysabove10(C::ClimGrid)
  @argcheck C[9] == "tas"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 10, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = buildarray(C, dataout, numYears)

  # Return climGrid type containing the indice
  return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits="days", latunits=C.latunits, lonunits=C.lonunits, variable="daysabove10", typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    daysabove10(data::AbstractArray, timevector::StepRange{Date,Base.Dates.Day})

Annual number with temperature >= 10 Celsius. This function returns a boolean vector. *true* if the data is higher or equal to 1 and *false* otherwise.
"""

function daysabove10(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 10, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function daysabove10(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 10, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function daysabove10(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 10, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

# """
#     customthresover(C::ClimGrid)
#
# customthresover, annual number of days over a specified threshold.
#
# Let TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:
#
#   TS[i,j] > thres.
# """
#
# function customthresover(C::ClimGrid, thres)
#   years    = Dates.year.(C.data[Axis{:time}][:])
#   numYears = unique(years)
#   dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
#   datain   = C.data.data
#
#   # Indice calculation
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t > thres, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
#   end
#
#   # Build output AxisArray
#   FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))
#
#   # Return climGrid type containing the indice
#   return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, variable = string("Days over ", thres, " ", C.dataunits), typeofvar = C.typeofvar, typeofcal = C.typeofcal)
# end
#
# """
#     customthresover(data::AbstractArray, time::StepRange{Date,Base.Dates.Day}, thres)
#
# customthresover, annual number of days over a specified threshold.
#
# Let TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:
#
#   TS[i,j] > thres.
# """
#
# function customthresover(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
#   years    = Dates.year.(timeV)
#   numYears = unique(years)
#   FD       = zeros(Int64, (length(numYears)))
#
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t > thres, +, view(FD, i:i), view(data, idx))
#   end
#   return FD
# end
#
# function customthresover(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
#   years    = Dates.year.(timeV)
#   numYears = unique(years)
#   FD       = zeros(Int64, (length(numYears), size(data, 2)))
#
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t > thres, +, view(FD, i:i, :), view(data, idx, :, :))
#   end
#   return FD
# end
#
# function customthresover(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
#   years    = Dates.year.(timeV)
#   numYears = unique(years)
#   FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))
#
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t > thres, +, view(FD, i:i, :, :), view(data, idx, :, :))
#   end
#   return FD
# end
#
# """
#     customthresunder(C::ClimGrid)
#
# customthresover, annual number of days under a specified threshold.
#
# Let TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:
#
#   TS[i,j] < thres.
# """
#
# function customthresunder(C::ClimGrid, thres)
#   years    = Dates.year.(C.data[Axis{:time}][:])
#   numYears = unique(years)
#   dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
#   datain   = C.data.data
#
#   # Indice calculation
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t < thres, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
#   end
#
#   # Build output AxisArray
#   FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))
#
#   # Return climGrid type containing the indice
#   return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, variable = string("Days under ", thres, " ", C.dataunits), typeofvar = C.typeofvar, typeofcal = C.typeofcal)
# end
#
#
# """
#     customthresunder(data::AbstractArray, time::StepRange{Date,Base.Dates.Day}, thres)
#
# customthresover, annual number of days under a specified threshold.
#
# Let TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:
#
#     TS[i,j] < thres.
# """
#
# function customthresunder(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
#   years    = Dates.year.(timeV)
#   numYears = unique(years)
#   FD       = zeros(Int64, (length(numYears)))
#
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t < thres, +, view(FD, i:i), view(data, idx))
#   end
#   return FD
# end
#
# function customthresunder(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
#   years    = Dates.year.(timeV)
#   numYears = unique(years)
#   FD       = zeros(Int64, (length(numYears), size(data, 2)))
#
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t < thres, +, view(FD, i:i, :), view(data, idx, :, :))
#   end
#   return FD
# end
#
# function customthresunder(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
#   years    = Dates.year.(timeV)
#   numYears = unique(years)
#   FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))
#
#   Threads.@threads for i in 1:length(numYears)
#     idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
#     Base.mapreducedim!(t -> t < thres, +, view(FD, i:i, :, :), view(data, idx, :, :))
#   end
#   return FD
# end
