"""
    prcp1(data::Array, timevector::StepRange{Date,Base.Dates.Day})

Annual number with preciptation over 1 mm. This function returns a boolean vector. *true* if the data is higher or equal to 1 and *false* otherwise.
"""

function prcp1(C::ClimGrid)
  @argcheck C[9] == "pr"
  years    = Dates.year(C.data[Axis{:time}][:])
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
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, var = "prcp1", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

function prcp1(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 1, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function prcp1(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 1, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function prcp1(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 1, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
  frostdays(data::Array, time::StepRange{Date,Base.Dates.Day})

FD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 Celsius.

Let TN(i,j) be daily minimum temperature on day i in year j. Count the number of days where:

  TN(i,j) < 0 Celsius.
"""

function frostdays(C::ClimGrid)
  @argcheck C[9] == "tasmin"
  years    = Dates.year(C.data[Axis{:time}][:])
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
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, var = "frostdays", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

function frostdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function frostdays(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function frostdays(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
  summerdays(data::Array, time::StepRange{Date,Base.Dates.Day})

SD, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25 degree Celsius.

Let TX(i,j) be daily maximum temperature on day i in year j. Count the number of days where:

  TX(i,j) >= 25 Celsius.
"""

function summerdays(C::ClimGrid)
  @argcheck C[9] == "tasmax"
  years    = Dates.year(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 25, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, var = "summerdays", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

function summerdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 25, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function summerdays(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 25, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function summerdays(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 25, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
  icingdays(data::Array, time::StepRange{Date,Base.Dates.Day})

ID, Number of summer days: Annual count of days when TX (daily maximum temperature) < 0 degree Celsius.

Let TX(i,j) be daily maximum temperature on day i in year j. Count the number of days where:

  TX(i,j) < 0 Celsius.
"""

function icingdays(C::ClimGrid)
  @argcheck C[9] == "tasmax"
  years    = Dates.year(C.data[Axis{:time}][:])
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
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, var = "icingdays", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

function icingdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
  return frostdays(data, timeV)
end

function icingdays(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  return frostdays(data, timeV)
end

function icingdays(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day})
  return frostdays(data, timeV)
end

"""
  tropicalnights(data::Array, time::StepRange{Date,Base.Dates.Day})

TropicalNights, Number of tropical nights: Annual count of days when TN (daily maximum temperature) > 20 degree Celsius.

Let TN(i,j) be daily minimum temperature on day i in year j. Count the number of days where:

  TN(i,j) > 20 Celsius.
"""

function tropicalnights(C::ClimGrid)
  @argcheck C[9] == "tasmin"
  years    = Dates.year(C.data[Axis{:time}][:])
  numYears = unique(years)
  dataout  = zeros(Int64, (length(numYears), size(C.data, 2), size(C.data, 3)))
  datain   = C.data.data

  # Indice calculation
  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 20, +, view(dataout, i:i, :, :), view(datain, idx, :,:))
  end

  # Build output AxisArray
  FD = AxisArray(dataout, Axis{:time}(Dates.Year.(numYears)), Axis{:lon}(C[1][Axis{:lon}][:]), Axis{:lat}(C[1][Axis{:lat}][:]))

  # Return climGrid type containing the indice
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, var = "tropicalnights", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

function tropicalnights(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 20, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function tropicalnights(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 20, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function tropicalnights(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 20, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
  customthresover(data::Array, time::StepRange{Date,Base.Dates.Day}, thres)

customthresover, annual number of days over a specified threshold.

Let TS(i,j) be a daily time serie value on day i in year j. Count the number of days where:

  TS(i,j) > thres.
"""

function customthresover(C::ClimGrid, thres)
  years    = Dates.year(C.data[Axis{:time}][:])
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
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, var = string("Days over ", thres, " ", C.dataunits), typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

function customthresover(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > thres, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function customthresover(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > thres, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function customthresover(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t > thres, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
  customthresunder(data::Array, time::StepRange{Date,Base.Dates.Day}, thres)

customthresover, annual number of days under a specified threshold.

Let TS(i,j) be a daily time serie value on day i in year j. Count the number of days where:

  TS(i,j) < thres.
"""

function customthresunder(C::ClimGrid, thres)
  years    = Dates.year(C.data[Axis{:time}][:])
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
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = "days", latunits = C.latunits, lonunits = C.lonunits, var = string("Days under ", thres, " ", C.dataunits), typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

function customthresunder(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < thres, +, view(FD, i:i), view(data, idx))
  end
  return FD
end

function customthresunder(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < thres, +, view(FD, i:i, :), view(data, idx, :, :))
  end
  return FD
end

function customthresunder(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day}, thres)
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < thres, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

"""
  annualmax(data::Array, time::StepRange{Date,Base.Dates.Day})

AM, Value of annual maximum of array data.

Let data(i,j) be daily time serie on day i in year j. Extract the highest value for year j.
"""

function annualmax(C::ClimGrid)
  years    = Dates.year(C.data[Axis{:time}][:])
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
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, var = "annualmax", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

function annualmax(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.maximum!(view(FD, i:i), view(data, idx))
  end
  return FD
end

function annualmax(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.maximum!(view(FD,i:i,:), view(data,idx,:))
  end
  return FD
end

function annualmax(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.maximum!(view(FD,i:i,:,:), view(data,idx,:,:))
  end
  return FD
end

function annualmax(data::Array{Float64, 3}, timeV::Array{Date,1})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.maximum!(view(FD,i:i,:,:), view(data,idx,:,:))
  end
  return FD
end

"""
  annualmin(data::Array, time::StepRange{Date,Base.Dates.Day})

AM, Value of annual minimum of array data.

Let data(i,j) be daily time serie on day i in year j. Extract the lowest value for year j.
"""

function annualmin(C::ClimGrid)
  years    = Dates.year(C.data[Axis{:time}][:])
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
  return ClimGrid(FD, model = C.model, experiment = C.experiment, run = C.run, filename = C.filename, dataunits = C.dataunits, latunits = C.latunits, lonunits = C.lonunits, var = "annualmin", typeofvar = C.typeofvar, typeofcal = C.typeofcal)
end

function annualmin(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.minimum!(view(FD, i:i), view(data, idx))
  end
  return FD
end

function annualmin(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.minimum!(view(FD,i:i,:), view(data,idx,:))
  end
  return FD
end

function annualmin(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Float64, (length(numYears), size(data, 2), size(data, 3)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.minimum!(view(FD,i:i,:,:), view(data,idx,:,:))
  end
  return FD
end
