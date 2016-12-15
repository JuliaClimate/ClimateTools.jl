"
  prcp1(data::Array, timevector::StepRange{Date,Base.Dates.Day})

Annual number with preciptation over 1 mm. This function returns a boolean vector. *true* if the data is higher or equal to 1 and *false* otherwise."


function prcp1(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 1, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

function prcp1(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t >= 1, +, view(FD, i:i, :, :), view(data, idx, :, :))
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

"
  frostdays(data::Array, time::StepRange{Date,Base.Dates.Day})

FD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 Celsius.

Let TN(i,j) be daily minimum temperature on day i in year j. Count the number of days where:

  TN(i,j) < 0 Celsius."

function frostdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(FD, i:i, :, :), view(data, idx, :, :))
  end
  return FD
end

function frostdays(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  years    = Dates.year(timeV)
  numYears = unique(years)
  FD       = zeros(Int64, (length(numYears), size(data, 2)))

  Threads.@threads for i in 1:length(numYears)
    idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
    Base.mapreducedim!(t -> t < 0, +, view(FD, i:i, :, :), view(data, idx, :, :))
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

"
  summerdays(TX::Array, time::StepRange{Date,Base.Dates.Day})

SD, Number of summer days: Annual count of days when TX (daily maximum temperature) > 25 degree Celsius.

Let TX(i,j) be daily maximum temperature on day i in year j. Count the number of days where:

  TX(i,j) >= 25 Celsius."

  function summerdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
    years    = Dates.year(timeV)
    numYears = unique(years)
    FD       = zeros(Int64, (length(numYears), size(data, 2)))

    Threads.@threads for i in 1:length(numYears)
      idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
      Base.mapreducedim!(t -> t >= 25, +, view(FD, i:i, :, :), view(data, idx, :, :))
    end
    return FD
  end

  function summerdays(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
    years    = Dates.year(timeV)
    numYears = unique(years)
    FD       = zeros(Int64, (length(numYears), size(data, 2)))

    Threads.@threads for i in 1:length(numYears)
      idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
      Base.mapreducedim!(t -> t >= 25, +, view(FD, i:i, :, :), view(data, idx, :, :))
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

"
  icingdays(TX::Array, time::StepRange{Date,Base.Dates.Day})

ID, Number of summer days: Annual count of days when TX (daily maximum temperature) < 0 degree Celsius.

Let TX(i,j) be daily maximum temperature on day i in year j. Count the number of days where:

  TX(i,j) < 0 Celsius."

  function icingdays(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
    return frostdays(data, timeV)
  end

  function icingdays(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
    return frostdays(data, timeV)
  end

  function icingdays(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day})
    return frostdays(data, timeV)
  end

"
  tropicalnights(TN::Array, time::StepRange{Date,Base.Dates.Day})

TropicalNights, Number of tropical nights: Annual count of days when TN (daily maximum temperature) > 20 degree Celsius.

Let TN(i,j) be daily minimum temperature on day i in year j. Count the number of days where:

  TN(i,j) > 20 Celsius."

  function tropicalnights(data::Array{Float64, 1}, timeV::StepRange{Date, Base.Dates.Day})
    years    = Dates.year(timeV)
    numYears = unique(years)
    FD       = zeros(Int64, (length(numYears), size(data, 2)))

    Threads.@threads for i in 1:length(numYears)
      idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
      Base.mapreducedim!(t -> t >= 20, +, view(FD, i:i, :, :), view(data, idx, :, :))
    end
    return FD
  end

  function tropicalnights(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
    years    = Dates.year(timeV)
    numYears = unique(years)
    FD       = zeros(Int64, (length(numYears), size(data, 2)))

    Threads.@threads for i in 1:length(numYears)
      idx = searchsortedfirst(years, numYears[i]):searchsortedlast(years, numYears[i])
      Base.mapreducedim!(t -> t >= 20, +, view(FD, i:i, :, :), view(data, idx, :, :))
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
