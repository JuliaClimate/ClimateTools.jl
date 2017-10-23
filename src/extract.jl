"""
    nc2julia(file::String, var::String; poly::Array)

Returns a ClimGrid type with the data in *file* of variable *var* inside the polygon *poly*. Metadata is built-in the ClimGrid type.

Inside the ClimgGrid type, the data is stored into an AxisArray data type, with time, longitude and latitude dimensions.
"""

function nc2julia(file::String, variable::String; poly = [])
  # TODO Finish polygon feature

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
  if !isempty(poly)
    # TODO extract index inside polygon
    msk = inpolyvec(lon, lat, poly)
  end

  # Get rectangular limits of the polygon


  # Get Data
  data = NetCDF.open(file, variable)
  if !isempty(poly)
    data = data[msk', :]
  end
  # Extract variable over a given region
  if variable == "pr"
    data = data[:,:,:] * 86400
    dataunits = "mm/day"
  elseif variable == "tasmax" || variable == "tasmin" || variable == "tas"
    data = data[:,:,:]
  end

  # # Convert to Float64 if Float32
  # if typeof(data[1]) == Float32
    #   data2 = Float64.(data)
  #   data = convert.(Array{Float64, 3}, data)
  #   # Make sure lat and lon are also Float64
  #   lon = convert(Array{Float64, 1}, lon)
  #   lat = convert(Array{Float64, 1}, lat)
  # end

  if dataunits == "K"
    data = data - 273.15
    dataunits = "Celsius"
  end

  # Permute dims --> make the longest dimension at position #1 (calculations are usually faster)
  if ndims(data) == 3
    data = permutedims(data, [3, 1, 2])
    # Convert data to AxisArray
    dataOut = AxisArray(data, Axis{:time}(timeV), Axis{:lon}(lon), Axis{:lat}(lat))
  elseif ndims(data) == 4 # this imply a 3D field (height component)
    data = permutedims(data, [4, 1, 2, 3])
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
    sumleapyear(dateV::)
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
This function return the polygon contained in shp (return type of Shapefile.jl package). It returns an Array{Float64,2}.

    shp2poly(shp::Shapefile.Handle{Shapefile.Polygon{Float64}})

    N.B. right now, it only works on the 1st polygon contained in shapefile shp. If your shapefile contains more than 1 polygon, it is suggested to use the function poly2array() through a for loop and build a vector of polygons.

"""
function shp2poly(shp::Shapefile.Handle{Shapefile.Polygon{Float64}})

    # Test if Shapefile Polygon is wrapping completely (closed Polygon)

    cond = shp.shapes[1].points[1].x == shp.shapes[1].points[end].x && shp.shapes[1].points[1].y == shp.shapes[1].points[end].y

    if cond
        P = Array{Float64}(2, length(shp.shapes[1].points))
    else # i.e. will need to wrap the return polygon P
        P = Array{Float64}(2, length(shp.shapes[1].points) + 1)
    end

  for i = 1:length(shp.shapes[1].points)
    P[1, i] = shp.shapes[1].points[i].x
    P[2, i] = shp.shapes[1].points[i].y
    if i == length(shp.shapes[1].points) && !cond
      P[1, i + 1] = shp.shapes[1].points[1].x
      P[2, i + 1] = shp.shapes[1].points[1].y
    end
  end
  return P
end

"""
This function return the polygon contained in shp (return type of Shapefile.jl package). It returns an Array{Float64,2}.

    shp2poly(shp::Shapefile.Polygon{Float64})

"""
function shp2poly(shp::Shapefile.Polygon{Float64})

    # Test if Shapefile Polygon is wrapping completely (closed Polygon)

    cond = shp.points[1].x == shp.points[end].x && shp.points[1].y == shp.points[end].y

    if cond
        P = Array{Float64}(2, length(shp.points))
    else # i.e. will need to wrap the return polygon P
        P = Array{Float64}(2, length(shp.points) + 1)
    end

  for i = 1:length(shp.points)
    P[1, i] = shp.points[i].x
    P[2, i] = shp.points[i].y
    if i == length(shp.points) && !cond
      P[1, i + 1] = shp.points[1].x
      P[2, i + 1] = shp.points[1].y
    end
  end
  return P
end

# function inpolygon{T<:Number}(x:: T, y:: T, vx:: Vector{T}, vy:: Vector{T})
#     @assert length(vx) == length(vy)
#     c = false
#     j = length(vx)
#     @inbounds for i=1:length(vx)
#         if (((vy[i] <= y && y < vy[j]) ||
#             (vy[j] <= y && y < vy[i])) &&
#             (x < (vx[j] - vx[i]) * (y - vy[i]) / (vy[j] - vy[i]) + vx[i]))
#             c = !c
#         end
#         j = i
#     end
#     return c
# end

# function fillpoly!{T}(M::Matrix{T}, px::Vector{Int}, py::Vector{Int}, value::T)
#     @assert length(px) == length(py)
#     left, right = minimum(px), maximum(px)
#     top, bottom = minimum(py), maximum(py)
#     @inbounds for x=left:right
#         ys = Set{Int64}()
#         j = length(px)
#         for i=1:length(px)
#             if (px[i] <= x && x <= px[j]) || (px[j] <= x && x <= px[i])
#                 # special case: adding the whole cut to ys
#                 if px[i] == px[j]
#                     push!(ys, py[i])
#                     push!(ys, py[j])
#                 else
#                     y = py[i] + (x - px[i]) / (px[j] - px[i]) * (py[j] - py[i])
#                     push!(ys, int(y))
#                 end
#             end
#             j = i
#         end
#         ys = sort([y for y in ys])
#         # if there's an odd number of intersection points, add one imeginary point
#         if length(ys) % 2 == 1
#             push!(ys, ys[end])
#         end
#         for i=1:2:length(ys)
#             M[ys[i]:ys[i+1], x] = value
#         end
#     end
#     return M
# end
#
# function poly2mask(px::Vector{Int}, py::Vector{Int}, m::Int, n::Int)
#     mask = zeros(Bool, m, n)
#     fillpoly!(mask, px, py, true)
# end
