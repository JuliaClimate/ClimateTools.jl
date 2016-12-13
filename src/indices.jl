"
  prcp1(data::Array, timevector::StepRange{Date,Base.Dates.Day})

Annual number with preciptation over 1 mm. This function returns a boolean vector. *true* if the data is higher or equal to 1 and *false* otherwise."


function prcp1(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  numYears = unique(Dates.year(timeV))
  PRCP1 = Array{Int64}(size(numYears, 1), size(data, 2))
  z = 1
  for iyear = numYears[1]:numYears[end]
    fgYear = findin(Dates.year(timeV), iyear)
    PRCP1[z, :] = sum([!isless(istep, 1) for istep in data[fgYear,:]], 1)
    z = z + 1
  end
  return PRCP1
  # return sum([isless(istep, 0) for istep in data], 1)
end

function prcp1(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day})
  numYears = unique(Dates.year(timeV))
  PRCP1 = Array{Int64}(size(numYears, 1), size(data, 2), size(data, 3))
  z = 1
  for iyear = numYears[1]:numYears[end]
    fgYear = findin(Dates.year(timeV), iyear)
    PRCP1[z, :, :] = sum([!isless(istep, 1) for istep in data[fgYear,:]], 1)
    z = z + 1
  end
  return PRCP1
  # return sum([isless(istep, 0) for istep in data], 1)
end
"
  frostdays(data::Array, time::StepRange{Date,Base.Dates.Day})

FD, Number of frost days: Annual count of days when TN (daily minimum temperature) < 0 Celsius.

Let TN(i,j) be daily minimum temperature on day i in year j. Count the number of days where:

  TN(i,j) < 0 Celsius."

function frostdays(data::Array{Float64, 2}, timeV::StepRange{Date, Base.Dates.Day})
  numYears = unique(Dates.year(timeV))
  FD = Array{Int64}(size(numYears, 1), size(data, 2))
  z = 1
  for iyear = numYears[1]:numYears[end]
    fgYear = findin(Dates.year(timeV), iyear)
    FD[z, :] = sum([isless(istep, 0) for istep in data[fgYear,:]], 1)
    z = z + 1
  end
  return FD
end

function frostdays(data::Array{Float64, 3}, timeV::StepRange{Date, Base.Dates.Day})
  numYears = unique(Dates.year(timeV))
  FD = Array{Int64}(size(numYears, 1), size(data, 2), size(data, 3))
  z = 1
  for iyear = numYears[1]:numYears[end]
    fgYear = findin(Dates.year(timeV), iyear)
    FD[z, :, :] = sum([isless(istep, 0) for istep in data[fgYear,:]], 1)
    z = z + 1
  end
  return FD
end
