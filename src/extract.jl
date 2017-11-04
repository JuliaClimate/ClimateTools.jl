"""
    nc2julia(file::String, var::String; poly::Array)

Returns a ClimGrid type with the data in *file* of variable *var* inside the polygon *poly*. Metadata is built-in the ClimGrid type.

Inside the ClimgGrid type, the data is stored into an AxisArray data type, with time, longitude and latitude dimensions.
"""

function nc2julia(file::String, variable::String; poly = Array{Float64}([]), data_units::String = "")

  # Get attributes for type "ClimGrid"
  ncI = NetCDF.ncinfo(file)
  experiment = NetCDF.ncgetatt(file, "global", "experiment_id")
  if !isa(experiment, String)
    experiment = "N/A"
  end
  runsim = NetCDF.ncgetatt(file, "global", "parent_experiment_rip")
  if !isa(runsim, String)
    runsim = "N/A"
  end
  model = NetCDF.ncgetatt(file, "global", "model_id")
  if !isa(model, String)
    model = "N/A"
  end
  dataunits = NetCDF.ncgetatt(file, variable, "units")
  latunits = NetCDF.ncgetatt(file, "lat", "units")
  lonunits = NetCDF.ncgetatt(file, "lon", "units")
  caltype = NetCDF.ncgetatt(file, "time", "calendar")

  # Construct time vector from info in netCDF file *str*
  timeV = buildtimevec(file)

  # Get index of grid points inside the polygon *poly*
  lat = NetCDF.ncread(file, "lat")
  lon = NetCDF.ncread(file, "lon")

  # Convert to positive values if lonunits are marked as "degrees_east" and they are all negative values
  if lonunits == "degrees_east" && sum(lon .< 0) == length(lon)
      lon += 360
  end

  # Get Data
  data = NetCDF.open(file, variable)
  if !isempty(poly)
    # TODO More consistent use of mask (the nc2julia function should return only the points inside the polygon. Other values should be set to NaN)
    msk = inpolyvec(lon, lat, poly)
    idlon, idlat = findn(.!isnan.(msk))
    minXgrid = minimum(idlon)
    maxXgrid = maximum(idlon)
    minYgrid = minimum(idlat)
    maxYgrid = maximum(idlat)
    if ndims(data) == 3
        data = data[minXgrid:maxXgrid, minYgrid:maxYgrid, :]
        # Permute dims
        data = permutedims(data, [3, 1, 2])
    elseif ndims(data) == 4
        data = data[minXgrid:maxXgrid, minYgrid:maxYgrid, :, :]
        # Permute dims
        data = permutedims(data, [4, 1, 2, 3])
    end

    #new mask (e.g. representing the region of the polygon)
    msk = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]
    data = applymask(data, msk)

    # Get lon-lat for such region
    lon = lon[minXgrid:maxXgrid]
    lat = lat[minYgrid:maxYgrid]

elseif isempty(poly) # no polygon clipping
    if ndims(data) == 3
        data = data[:, :, :]
        # Permute dims
        data = permutedims(data, [3, 1, 2])

    elseif ndims(data) == 4
        data = data[:, :, :, :]
        # Permute dims
        data = permutedims(data, [4, 1, 2, 3])
    end
  end

  # Convert units of optional argument data_units is provided
  if data_units == "Celsius" && (variable == "tas" || variable == "tasmax" || variable == "tasmin") && dataunits == "K"
    data = data - 273.15
    dataunits = "Celsius"
  end

  if data_units == "mm" && variable == "pr" && (dataunits == "kg m-2 s-1" || dataunits == "mm s-1")

    rez = timeresolution(file)
    factor = pr_timefactor(rez)
    data = data * factor
    if rez != "N/A"
        dataunits = string("mm/",rez)
    else
        dataunits = "mm"
    end
  end

  # Create AxisArray from variable "data"
  if ndims(data) == 3
    # Convert data to AxisArray
    dataOut = AxisArray(data, Axis{:time}(timeV), Axis{:lon}(lon), Axis{:lat}(lat))
  elseif ndims(data) == 4 # this imply a 3D field (height component)
    # Get level vector
    plev = NetCDF.ncread(file, "plev")
    # Convert data to AxisArray
    dataOut = AxisArray(data, Axis{:time}(timeV), Axis{:lon}(lon), Axis{:lat}(lat), Axis{:plev}(plev))
  else
    throw(error("nc2julia takes only 3D and 4D variables for the moment"))
  end

  return ClimGrid(dataOut, model = model, experiment = experiment, run = runsim, filename = file, dataunits = dataunits, latunits = latunits, lonunits = lonunits, variable = variable, typeofvar = variable, typeofcal = caltype)


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
    timeRaw = floor.(NetCDF.ncread(str, "time"))
    leapDaysPer = sumleapyear(initDate, timeRaw[1])
    leapDaysPer2 = sumleapyear(initDate, timeRaw[end])
    startDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[1]))) + Base.Dates.Day(leapDaysPer)
    endDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[end]))) + Base.Dates.Day(leapDaysPer2)

    dateTmp = Date(startDate):Date(endDate)
    # REMOVE leap year
    idx = (Dates.month.(dateTmp) .== 2) .&  (Dates.day.(dateTmp) .== 29)
    if length(idx) !== 1
      dateTmp = dateTmp[.!idx]
    else
      dateTmp = Array{Date}(dateTmp)
    end
elseif calType == "gregorian" || calType == "standard"
    timeRaw = floor.(NetCDF.ncread(str, "time"))
    leapDaysPer = sumleapyear(initDate, timeRaw[1])
    leapDaysPer2 = sumleapyear(initDate, timeRaw[end])
    startDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[1])))
    endDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[end])))
    dateTmp = Date(startDate):Date(endDate)
  end
  # output date vector
  dateTmp = convert(Array{Date, 1}, dateTmp)
  return dateTmp
end

"""
Number of leap years in date vector

    sumleapyear(dates::StepRange{Date,Base.Dates.Day})

    sumleapyear(initDate::Date, timeRaw)
"""
function sumleapyear(initDate::Date, timeRaw)

  out = 0::Int
  endDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[1])))
  years = unique(Dates.year.(initDate:endDate))
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
  years = unique(Dates.year.(dates))
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
This function return the polygons contained in shp.shapes[i] (return type of Shapefile.jl package). It returns the x and y coordinates vectors.

    shapefile_coords(poly::Shapefile.Polygon)

"""
function shapefile_coords(poly::Shapefile.Polygon)
    start_indices = poly.parts+1
    end_indices = vcat(poly.parts[2:end], length(poly.points))
    x, y = zeros(0), zeros(0)
    for (si,ei) in zip(start_indices, end_indices)
        push!(x, NaN)
        push!(y, NaN)
        for pt in poly.points[si:ei]
            push!(x, pt.x)
            push!(y, pt.y)
        end
    end
    x, y
end

"""
This function return the time resolution of the netCDF file "str"

    function timeresolution(str:String)

"""

function timeresolution(str::String)

    timevec = (NetCDF.ncread(str, "time"))
    if length(timevec) > 1
        timediff = diff(timevec)[1]
        if timediff == 1. || timediff == 1
            return "24h"
        elseif timediff == 0.5
            return "12h"
        elseif timediff == 0.25
            return "6h"
        elseif timediff == 0.125
            return "3h"
        end
    else
        return "N/A"
    end
end

"""
This function return the time factor that should be applied to ptrecipitation to get accumulation for resolution "rez"

    function pr_timefactor(rez:String)

"""

function pr_timefactor(rez::String)

    if rez == "24h"
        return 86400.
    elseif rez == "12h"
        return 43200.
    elseif rez == "6h"
        return 21600.
    elseif rez == "3h"
        return 10800.
    elseif rez == "N/A"
        return 1.
    end

end
