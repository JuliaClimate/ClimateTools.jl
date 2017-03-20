"""
    nc2julia(file::String, var::String; poly::Array{Float64})

Returns a ClimGrid type with the data in *file* of variable *var* inside the polygon *poly*. Metadata is built-in the ClimGrid type.

Inside the ClimgGrid type, the data is stored into an AxisArray data type, with time, longitude and latitude dimensions.
"""

function nc2julia(file::String, var::String; poly::Array{Float64} = [])
  # TODO Finish polygon feature

  # Get attributes for type "ClimGrid"
  ncI = NetCDF.ncinfo(file)
  experiment = NetCDF.ncgetatt(file, "global", "experiment_id")
  if !isa(experiment, String)
    experiment = ""
  end
  runsim = NetCDF.ncgetatt(file, "global", "parent_experiment_rip")
  if !isa(runsim, String)
    runsim = ""
  end
  model = NetCDF.ncgetatt(file, "global", "model_id")
  if !isa(model, String)
    model = ""
  end
  dataunits = NetCDF.ncgetatt(file, var, "units")
  latunits = NetCDF.ncgetatt(file, "lat", "units")
  lonunits = NetCDF.ncgetatt(file, "lon", "units")
  caltype = NetCDF.ncgetatt(file, "time", "calendar")

  # Construct time vector from info in netCDF file *str*
  timeV = buildtimevec(file)

  # Get index of grid points inside the polygon *poly*
  lat = NetCDF.ncread(file, "lat")
  lon = NetCDF.ncread(file, "lon")
  if !isempty(poly)
    # TODO extract index inside polygon
    # polyV = inpolyV(lat, lon, poly)
  end


  # Get Data
  data = NetCDF.open(file, var)
  if !isempty(poly)
    # TODO subset of whole data
  end
  if var == "pr"
    data = data[:,:,:] * 86400
    dataunits = "mm/day"
  elseif var == "tasmax" || var == "tasmin" || var == "tas"
    data = data[:,:,:]
  end

  # Convert to Float64 if Float32
  if typeof(data[1]) == Float32
    data = convert(Array{Float64}, data)
  end

  if dataunits == "K"
    data = data - 273.15
    dataunits = "Celsius"
  end

  # Permute dims --> make the longest dimension at position #1 (calculations are usually faster)
  data = permutedims(data, [3, 1, 2])
  # dataOut = AxisArray(data, :time, :lon, :lat)
  dataOut = AxisArray(data, Axis{:time}(timeV), Axis{:lon}(lon), Axis{:lat}(lat))#, :lon,:lat)

  return ClimGrid(dataOut, model = model, experiment = experiment, run = runsim, filename = file, dataunits = dataunits, latunits = latunits, lonunits = lonunits, var = var, typeofvar = var, typeofcal = caltype)


end

"""
    buildtimevec(str::String)

Construct the time vector from the netCDF file str

"""

function buildtimevec(str::String)

  # Time units
  units = NetCDF.ncgetatt(str, "time", "units") # get starting date
  m = match(r"(\d+)[-.\/](\d+)[-.\/](\d+)", units, 1) # match a date from string
  daysfrom = m.match # get only the date ()"yyyy-mm-dd" format)
  initDate = Date(daysfrom, "yyyy-mm-dd")

  # Calendar type
  calType = NetCDF.ncgetatt(str, "time", "calendar")
  if calType == "noleap" || calType == "365_day"
    nDays = 365
    # get time of netCDF file *str*
    timeRaw = floor(NetCDF.ncread(str, "time"))
    leapDaysPer = sumleapyear(initDate, timeRaw[1])
    leapDaysPer2 = sumleapyear(initDate, timeRaw[end])
    startDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[1]))) + Base.Dates.Day(leapDaysPer)
    endDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[end]))) + Base.Dates.Day(leapDaysPer2)

    dateTmp = Date(startDate):Date(endDate)
    # REMOVE leap year
    idx = Dates.monthday(dateTmp) .== (2,29)
    dateTmp = dateTmp[!idx]
  elseif calType == "gregorian"
    timeRaw = floor(NetCDF.ncread(str, "time"))
    leapDaysPer = sumleapyear(initDate, timeRaw[1])
    leapDaysPer2 = sumleapyear(initDate, timeRaw[end])
    startDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[1])))
    endDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[end])))
    dateTmp = Date(startDate):Date(endDate)
  end
  # output date vector
  return dateTmp
end

"""
    sumleapyear(dateV::)
"""
function sumleapyear(initDate::Date, timeRaw)

  out = 0::Int
  endDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[1])))
  years = unique(Dates.year(initDate:endDate))
  # Sum over time vector
  for idx = 1:length(years)
    if Dates.isleapyear(years[idx])
      out += 1
    end
    # out[idx] = Dates.isleapyear(test[idx])
  end

  return out

end

function sumleapyear(dates::StepRange{Date,Base.Dates.Day})

  out = 0::Int
  endDate = dates[end]#initDate + Base.Dates.Day(convert(Int64,round(timeRaw[1])))
  years = unique(Dates.year(dates))
  # Sum over time vector
  for idx = 1:length(years)
    if Dates.isleapyear(years[idx])
      out += 1
    end
    # out[idx] = Dates.isleapyear(test[idx])
  end

  return out

end

"""
    inpolyV(lat, lon, poly::Array{Float64})

"""

# function inpolyV(lat, lon, shp::String)
#   # TODO finish this feature
#   # Build polygon from shapefile *shp*
#   poly = shpextract(shp)
#   if !isempty(lat)
#     # dummy call for linting
#   elseif !isempty(lon)
#     # dummy call for linting
#   end
#
# end
#
# """
#     shpextract(shp::String)
# """
#
# function shpextract(shp::String)
#   shp = Shapefile.read(shp, Shapefile.Handle)
#   P = Array(Float64, 2, length(shp.shapes[2].points) + 1)
#   # CREATE VECTOR
#   for i = 1:length(shp.shapes[2].points)
#     P[1, i] = shp.shapes[2].points[i].x
#     P[2, i] = shp.shapes[2].points[i].y
#     if i == length(shp.shapes[2].points)
#       P[1, i + 1] = shp.shapes[2].points[1].x
#       P[2, i + 1] = shp.shapes[2].points[1].y
#     end
#   end
#
#
# end
