"""
    nc2julia(file::String, var::String, poly::Vector)

Returns a ClimGrid type with the data in *file* of variable *var* inside the polygon *poly*. Metadata is built-in the ClimGrid type.

Inside the ClimgGrid type, the data is stored into an AxisArray data type, with time, longitude and latitude dimensions.
"""

function nc2julia(file::String, var::String, poly::Array{Float64})

  # Get attributes for type "ClimGrid"
  ncI = NetCDF.ncinfo(file)
  experiment = NetCDF.ncgetatt(file, "global", "experiment_id")
  run = NetCDF.ncgetatt(file, "global", "parent_experiment_rip")
  model = NetCDF.ncgetatt(file, "global", "model_id")
  dataunits = NetCDF.ncgetatt(file, var, "units")
  latunits = NetCDF.ncgetatt(file, "lat", "units")
  lonunits = NetCDF.ncgetatt(file, "lon", "units")
  caltype = NetCDF.ncgetatt(file, "time", "calendar")

  # Construct time vector from info in netCDF file *str*
  timeV = buildtimevec(file)

  # Get index of grid points inside the polygon *poly*
  lat = NetCDF.ncread(file, "lat")
  lon = NetCDF.ncread(file, "lon")
  # polyV = inpolyV(str, poly) #to-do!!

  # Get Data
  data = NetCDF.open(file, var)
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

  return ClimGrid(dataOut, model, experiment, run, file, dataunits, latunits, lonunits, var)


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
  if calType == "noleap"
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
  else
    # TO-DO!
    error("TO-DO!")
  end


  # output date vector
  return dateTmp


end

"""
    sumleapyear(dateV::)
"""
function sumleapyear(initDate::Date, timeRaw::Real)

  out = 0::Int64
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

"""
    inpolyV(str::String, poly::Array{Float64})

"""

function inpolyV(str::String, shp::String)
  # Build polygon from shapefile *shp*
  poly = shpextract(shp)

end

"""
    shpextract(shp::String)
"""

function shpextract(shp::String)
  shp = Shapefile.read(shp, Shapefile.Handle)
  P = Array(Float64, 2, length(shp.shapes[2].points) + 1)
  # CREATE VECTOR
  for i = 1:length(shp.shapes[2].points)
    P[1, i] = shp.shapes[2].points[i].x
    P[2, i] = shp.shapes[2].points[i].y
    if i == length(shp.shapes[2].points)
      P[1, i + 1] = shp.shapes[2].points[1].x
      P[2, i + 1] = shp.shapes[2].points[1].y
    end
  end


end
