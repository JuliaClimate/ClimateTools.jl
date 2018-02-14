"""
    prcp1(C::ClimGrid)

Annual number with preciptation >= 1 mm. This function returns a ClimGrid.
"""

function prcp1(C::ClimGrid)
  @argcheck C[9] == "pr"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 1, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, variable = "prcp1", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

"""
    prcp1(data::AbstractArray, timevector::StepRange{Date,Base.Dates.Day})

Annual number with preciptation >= 1 mm. This function returns a boolean vector. *true* if the data is higher or equal to 1 and *false* otherwise.
"""

function prcp1(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 1, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function prcp1(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 1, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function prcp1(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 1, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
    frostdays(C::ClimGrid)

FD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 Celsius.

Let TN[i,j] be daily minimum temperature on day i in year j. Count the number of days where:

  TN[i,j] < 0 Celsius.
"""

function frostdays(C::ClimGrid)
  @argcheck C[9] == "tasmin"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, variable = "frostdays", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

"""
    frostdays(data::AbstractArray, time::StepRange{Date,Base.Dates.Day})

FD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 Celsius.

Let TN[i,j] be daily minimum temperature on day i in year j. Count the number of days where:

  TN[i,j] < 0 Celsius.
"""

function frostdays(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function frostdays(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function frostdays(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
    summerdays(C::ClimGrid)

SD, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25 degree Celsius.

Let TX[i,j] be daily maximum temperature on day i in year j. Count the number of days where:

  TX[i,j] > 25 Celsius.
"""

function summerdays(C::ClimGrid)
  @argcheck C[9] == "tasmax"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > 25, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, variable = "summerdays", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end


"""
    summerdays(data::AbstractArray, time::StepRange{Date,Base.Dates.Day})

SD, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25 degree Celsius.

Let TX[i,j] be daily maximum temperature on day i in year j. Count the number of days where:

  TX[i,j] > 25 Celsius.
"""

function summerdays(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > 25, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function summerdays(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > 25, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function summerdays(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > 25, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
    icingdays(C::ClimGrid)

ID, Number of summer days: Annual count of days when TX (daily maximum temperature) < 0 degree Celsius.

Let TX[i,j] be daily maximum temperature on day i in year j. Count the number of days where:

  TX[i,j] < 0 Celsius.
"""

function icingdays(C::ClimGrid)
  @argcheck C[9] == "tasmax"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, variable = "icingdays", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end


"""
    icingdays(data::AbstractArray, time::StepRange{Date,Base.Dates.Day})

ID, Number of summer days: Annual count of days when TX (daily maximum temperature) < 0 degree Celsius.

Let TX[i,j] be daily maximum temperature on day i in year j. Count the number of days where:

  TX[i,j] < 0 Celsius.
"""

function icingdays(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
  return frostdays(data, timeV)
end

function icingdays(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
  return frostdays(data, timeV)
end

function icingdays(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
  return frostdays(data, timeV)
end

"""
    tropicalnights(C::ClimGrid)

TropicalNights, Number of tropical nights: Annual count of days when TN (daily maximum temperature) > 20 degree Celsius.

Let TN[i,j] be daily minimum temperature on day i in year j. Count the number of days where:

  TN[i,j] > 20 Celsius.
"""

function tropicalnights(C::ClimGrid)
  @argcheck C[9] == "tasmin"
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > 20, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, variable = "tropicalnights", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

"""
    tropicalnights(data::AbstractArray, time::StepRange{Date,Base.Dates.Day})

TropicalNights, Number of tropical nights: Annual count of days when TN (daily maximum temperature) > 20 degree Celsius.

Let TN[i,j] be daily minimum temperature on day i in year j. Count the number of days where:

  TN[i,j] > 20 Celsius.
"""

function tropicalnights(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > 20, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function tropicalnights(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > 20, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function tropicalnights(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > 20, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
    customthresover(C::ClimGrid)

customthresover, annual number of days over a specified threshold.

Let TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:

  TS[i,j] > thres.
"""

function customthresover(C::ClimGrid, thres)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > thres, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, variable = string("Days over ", thres, " ", C.dataunits), typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

"""
    customthresover(data::AbstractArray, time::StepRange{Date,Base.Dates.Day}, thres)

customthresover, annual number of days over a specified threshold.

Let TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:

  TS[i,j] > thres.
"""

function customthresover(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > thres, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function customthresover(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > thres, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function customthresover(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > thres, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
    customthresunder(C::ClimGrid)

customthresover, annual number of days under a specified threshold.

Let TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:

  TS[i,j] < thres.
"""

function customthresunder(C::ClimGrid, thres)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < thres, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, variable = string("Days under ", thres, " ", C.dataunits), typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end


"""
    customthresunder(data::AbstractArray, time::StepRange{Date,Base.Dates.Day}, thres)

customthresover, annual number of days under a specified threshold.

Let TS[i,j] be a daily time serie value on day i in year j. Count the number of days where:

    TS[i,j] < thres.
"""

function customthresunder(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < thres, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function customthresunder(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < thres, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function customthresunder(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < thres, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
    annualmax(C::ClimGrid)

Annual maximum of array data.

Let data[i,j] be daily time serie on day i in year j. Extract the highest value for year j.
"""

function annualmax(C::ClimGrid)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Float64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.maximum!(view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, variable = "annualmax", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

"""
    annualmax(data::AbstractArray, time::StepRange{Date,Base.Dates.Day})

Annual maximum of array data.

Let data[i,j] be daily time serie on day i in year j. Extract the highest value for year j.
"""

function annualmax(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.maximum!(view(FD, i:i), view(data, idx))
  end
  return FD
end

function annualmax(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.maximum!(view(FD,i:i,:), view(data,idx,:))
  end
  return FD
end

function annualmax(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.maximum!(view(FD,i:i,:,:), view(data,idx,:,:))
  end
  return FD
end

function annualmax(data::AbstractArray{N, 3} where N, timeV::Array{Date,1})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.maximum!(view(FD,i:i,:,:), view(data,idx,:,:))
  end
  return FD
end

"""
    annualmin(C::ClimGrid)

Annual minimum of array data.

Let data[i,j] be daily time serie on day i in year j. Extract the lowest value for year j.
"""

function annualmin(C::ClimGrid)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Float64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.minimum!(view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, variable = "annualmin", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

"""
    annualmin(data::AbstractArray, time::StepRange{Date,Base.Dates.Day})

Annual minimum of array data.

Let data[i,j] be daily time serie on day i in year j. Extract the lowest value for year j.
"""

function annualmin(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.minimum!(view(FD, i:i), view(data, idx))
  end
  return FD
end

function annualmin(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.minimum!(view(FD,i:i,:), view(data,idx,:))
  end
  return FD
end

function annualmin(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.minimum!(view(FD,i:i,:,:), view(data,idx,:,:))
  end
  return FD
end

"""
    annualmean(C::ClimGrid)

Annual mean of array data.

Let data[i,j] be daily time serie on day i in year j. Calculate the mean value for year j.
"""
function annualmean(C::ClimGrid)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Float64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mean!(view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, variable = "annualmean", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

"""
    annualmean(data::AbstractArray, time::StepRange{Date,Base.Dates.Day})

Annual mean of array data.

Let data[i,j] be daily time serie on day i in year j. Calculate the mean value for year j.
"""

function annualmean(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mean!(view(FD, i:i), view(data, idx))
  end
  return FD
end

function annualmean(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mean!(view(FD,i:i,:), view(data,idx,:))
  end
  return FD
end

function annualmean(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mean!(view(FD,i:i,:,:), view(data,idx,:,:))
  end
  return FD
end

"""
    annualsum(C::ClimGrid)

Annual sum of array data.

Let data[i,j] be daily time serie on day i in year j. Sums daily values for year j.
"""
function annualsum(C::ClimGrid)
  years    = Dates.year.(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Float64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.sum!(view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, variable = "annualsum", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

"""
    annualsum(data::AbstractArray, time::StepRange{Date,Base.Dates.Day})

Value of annual sum of array data.

Let data[i,j] be daily time serie on day i in year j. Sums daily values for year j.
"""

function annualsum(data::AbstractArray{N, 1} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.sum!(view(FD, i:i), view(data, idx))
  end
  return FD
end

function annualsum(data::AbstractArray{N, 2} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.sum!(view(FD,i:i,:), view(data,idx,:))
  end
  return FD
end

function annualsum(data::AbstractArray{N, 3} where N, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year.(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.sum!(view(FD,i:i,:,:), view(data,idx,:,:))
  end
  return FD
end
