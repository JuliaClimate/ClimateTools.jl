"""
    nc2julia(file::String, variable::String; poly = Array{Float64}([]), start_date::Date, end_date::Date, data_units::String = "")

Returns a ClimGrid type with the data in **file** of variable **var** inside the polygon **poly**. Metadata is built-in the ClimGrid type, from the netCDF attributes.

Inside the ClimgGrid type, the data is stored into an AxisArray data type, with time, longitude/x and latitude/y dimensions.

The polygon provided should be in the -180, +180 longitude format. If the polygon crosses the International Date Line, the polygon should be splitted in multiple parts (i.e. multi-polygons).

Options for data_units are for precipitation : "mm", which converts the usual "kg m-2 s-1" unit found in netCDF files. For temperature : "Celsius", which converts the usual "Kelvin" unit.

**Note:** nc2julia is based on CF conventions (http://cfconventions.org/). If you are unable to read the netCDF file with nc2julia, the user will need to read it with low-level functions available in the NetCDF.jl package (https://github.com/JuliaGeo/NetCDF.jl).
"""

function nc2julia(file::String, variable::String; poly = ([]), start_date::Date = Date(-4000), end_date::Date = Date(-4000), data_units::String = "")


  # TODO this file is a complete mess, but it works. Clean it up!
  # TODO create another function for orography extraction (2D dataset)
  # Get attributes for type "ClimGrid" using NCDataset.jl
  ncI = NetCDF.ncinfo(file);
  attribs = NetCDF.ncinfo(file).gatts;
  ds = NCDatasets.Dataset(file)
  # ncI = NetCDF.ncinfo(file);
  attribs_dataset = ds.attrib

  # The following attributes should be set for netCDF files that follows CF conventions.
  # project, institute, model, experiment, frequency

  project = project_id(attribs_dataset)
  institute = institute_id(attribs_dataset)
  model = model_id(attribs_dataset)
  experiment = experiment_id(attribs_dataset)
  frequency = frequency_var(attribs_dataset)
  runsim = runsim_id(attribs_dataset)

  dataunits = ds[variable].attrib["units"]
  latunits = ds["lat"].attrib["units"]
  lonunits = ds["lon"].attrib["units"]
  caltype = ds["time"].attrib["calendar"]

  # Get dimensions names
  latname = getdim_lat(ds)
  lonname = getdim_lon(ds)

  # Create dict with latname and lonname
  dimension_dict = Dict(["lon" => lonname, "lat" => latname])

  lat_raw = NetCDF.ncread(file, latname)
  lon_raw = NetCDF.ncread(file, lonname)

  # Get variable attributes
  varattrib = Dict(ds[variable].attrib)

  if latname != "lat" && lonname != "lon" # means we don't have a "regular" grid
      latgrid = NetCDF.ncread(file, "lat")#ds["lat"][:, :]
      longrid = NetCDF.ncread(file, "lon")#ds["lon"][:, :]
      if @isdefined varattrib["grid_mapping"]
          map_dim = varattrib["grid_mapping"]
          map_attrib = Dict(ds[map_dim].attrib)
          map_attrib["grid_mapping"] = map_dim
      else
          error("File an issue on https://github.com/Balinus/ClimateTools.jl/issues to get the grid supported")
      end
  else # if no grid provided, create one
      # longrid, latgrid = ndgrid(lon_raw, lat_raw)
      # longrid, latgrid = meshgrid(lon_raw, lat_raw)
      longrid, latgrid = ndgrid(lon_raw, lat_raw)
      map_attrib = Dict(["grid_mapping_name" => "Regular_longitude_latitude"])
  end

  # =====================
  # TIME

  # Construct time vector from info in netCDF file *str*
  # TODO add time frequency to buildtimevec to comply with monthly values
  timeV = buildtimevec(file)
  if frequency == "mon"
      timeV = corr_timevec(timeV, frequency)
  end

  if start_date !== Date(-4000)
      @argcheck start_date <= end_date
      @argcheck start_date >= timeV[1]
      @argcheck end_date <= timeV[end]
      if frequency != "mon"
          idxtimebeg = find(timeV .== start_date)[1]
          idxtimeend = find(timeV .== end_date)[1]
      elseif frequency == "mon"
          idxtimebeg = find(timeV .== start_date)[1]
          idxtimeend = find(timeV .== Date(Dates.year(end_date), Dates.month(end_date), Dates.day(1)))[1]

      end
  else
      idxtimebeg = 1
      idxtimeend = length(timeV)
  end

  timeV = timeV[idxtimebeg:idxtimeend]

  # Correct longitude from 0, 360 => -180, 180
  # from 0, 360 to -180, 180
  # dataout, lonout = basemap[:shiftgrid](180.0, datain', longrid[:, 1], start=false)
  # from -180, 180 to 0, 360
  # dataorig, lonorig = basemap[:shiftgrid](0.0, dataout, lonout, start=true)

  rotatedgrid = false
  if sum(longrid .> 180) >= 1
      rotatedgrid = true

      for ilon in eachindex(longrid)
          if longrid[ilon] >= 180
              longrid[ilon] = longrid[ilon] - 360
          end
      end

      for ilon in eachindex(lon_raw)
          if lon_raw[ilon] >= 180
              lon_raw[ilon] = lon_raw[ilon] - 360
          end
      end

      ieast = longrid .>= 0
      iwest = longrid .< 0

      grideast = reshape(longrid[ieast], :, size(longrid, 2))
      gridwest = reshape(longrid[iwest], :, size(longrid, 2))

      longrid_flip = vcat(gridwest, grideast)

      ieast_vec = lon_raw .>= 0
      iwest_vec = lon_raw .< 0

      lon_raw_flip = [lon_raw[iwest_vec];lon_raw[ieast_vec]]


  else
      longrid_flip = longrid


  end



  # ===================
  # GET DATA
  # data = ds[variable]
  data = NetCDF.open(file, variable)
  if !isempty(poly)

    minpoly = minimum(poly[1, .!isnan.(poly[1, :])])
    maxpoly = maximum(poly[1, .!isnan.(poly[1, :])])

    meridian = false
    if sign(minpoly)*sign(maxpoly) == sign(-1.0) # polygon crosses the meridian
       meridian = true
    end

    msk = inpolygrid(longrid_flip, latgrid, poly) # mask of inpoly

    if sum(isnan.(msk)) == length(msk) # no grid point insode polygon
        throw(error("No grid points found inside the provided polygon"))
    end

    if rotatedgrid
        # Regrid msk if the grid have been rotated to get proper index for data extraction
        idxeast = longrid_flip .>= 0
        idxwest = longrid_flip .< 0
        mskeast = reshape(msk[idxeast], :, size(msk, 2))
        mskwest = reshape(msk[idxwest], :, size(msk, 2))

        msk = vcat(mskeast, mskwest)

    end

    idlon, idlat = findn(.!isnan.(msk))
    minXgrid = minimum(idlon)
    maxXgrid = maximum(idlon)
    minYgrid = minimum(idlat)
    maxYgrid = maximum(idlat)
    if ndims(data) == 3
        data = data[minXgrid:maxXgrid, minYgrid:maxYgrid, idxtimebeg:idxtimeend]
        # Permute dims
        # permutedims!(data, data, [3, 1, 2])
        data = permutedims(data, [3, 1, 2])
    elseif ndims(data) == 4
        data = data[minXgrid:maxXgrid, minYgrid:maxYgrid, :, idxtimebeg:idxtimeend]
        # Permute dims
        data = permutedims(data, [4, 1, 2, 3])
    end

    #new mask (e.g. representing the region of the polygon)
    msk = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]
    data = applymask(data, msk)

    # Get lon_raw and lat_raw for such region
    lon_raw = lon_raw[minXgrid:maxXgrid]
    lat_raw = lat_raw[minYgrid:maxYgrid]

    if map_attrib["grid_mapping_name"] == "Regular_longitude_latitude"

        ieastvec = lon_raw .>= 0
        iwestvec = lon_raw .< 0

        lon_raw = [lon_raw[iwestvec]; lon_raw[ieastvec]]
    end

    # Idem for longrid and latgrid
    if meridian
        longrid = vcat(grideast, gridwest)
        longrid = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
        latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]

        # Re-translate to -180, 180 if 0, 360

        if rotatedgrid

            idxeast = longrid .>= 0
            idxwest = longrid .< 0
            grideastsub = reshape(longrid[idxeast], :, size(longrid, 2))
            gridwestsub = reshape(longrid[idxwest], :, size(longrid, 2))

            longrid_flip = vcat(gridwestsub, grideastsub)

            data = permute_west_east(data, idxwest, idxeast)
        end

    else

        if rotatedgrid # flip in original grid
            longrid = vcat(grideast, gridwest)


        end

        longrid = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
        latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]

        # reflip in destination grid
        # if rotatedgrid # flip in original grid
        #
        # end


        if rotatedgrid
            # if sum(longrid .> 180) >= 1
                idxeast = longrid .>= 0
                idxwest = longrid .< 0
                grideastsub = reshape(longrid[idxeast], :, size(longrid, 2))
                gridwestsub = reshape(longrid[idxwest], :, size(longrid, 2))

                longrid_flip = vcat(gridwestsub, grideastsub)

                ieast_vec = lon_raw .>= 0
                iwest_vec = lon_raw .< 0

                lon_raw_flip = [lon_raw[ieast_vec];lon_raw[iwest_vec]]

                data = permute_west_east(data, idxwest, idxeast)
                # longrid = longrid .- 180.0
                # else
                # longrid = longrid + 180.0
            # else
                # longrid_flip = longrid
            # end


        end

    end

  elseif isempty(poly) # no polygon clipping
    if ndims(data) == 3
        data = data[:, :, idxtimebeg:idxtimeend]
        # # Permute dims (climate indices calculations are quicker, but extraction is longer)
        data = permutedims(data, [3, 1, 2])

    elseif ndims(data) == 4
        data = data[:, :, :, idxtimebeg:idxtimeend]
        # # Permute dims
        data = permutedims(data, [4, 1, 2, 3])
    end
    msk = Array{Float64}(ones((size(data, 2), size(data, 3))))

    if rotatedgrid
        # Flip data
        idxeast = longrid .>= 0
        idxwest = longrid .< 0
        data = permute_west_east(data, idxwest, idxeast)
    end

  end

  # # Replace fillvalues with NaN
  fillval = NetCDF.ncgetatt(file, variable, "_FillValue")
  data[data .== fillval] = NaN

  if rotatedgrid
      longrid = longrid_flip
      lon_raw = lon_raw_flip
  end

  # Convert units of optional argument data_units is provided
  if data_units == "Celsius" && (variable == "tas" || variable == "tasmax" || variable == "tasmin") && dataunits == "K"
    data = data - 273.15
    dataunits = "Celsius"
  end

  if data_units == "mm" && variable == "pr" && (dataunits == "kg m-2 s-1" || dataunits == "mm s-1")

    rez = timeresolution(NetCDF.ncread(file, "time"))
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
    dataOut = AxisArray(data, Axis{:time}(timeV), Axis{Symbol(lonname)}(lon_raw), Axis{Symbol(latname)}(lat_raw))
  elseif ndims(data) == 4 # this imply a 3D field (height component)
    # Get level vector
    plev = ds["plev"][:]#NetCDF.ncread(file, "plev")
    # Convert data to AxisArray
    dataOut = AxisArray(data, Axis{:time}(timeV), Axis{Symbol(lonname)}(lon_raw), Axis{Symbol(latname)}(lat_raw), Axis{:plev}(plev))
  else
    throw(error("nc2julia takes only 3D and 4D variables for the moment"))
  end

  close(ds)
  NetCDF.ncclose(file)

  return ClimGrid(dataOut, longrid=longrid, latgrid=latgrid, msk=msk, grid_mapping=map_attrib, dimension_dict=dimension_dict, model=model, frequency=frequency, experiment=experiment, run=runsim, project=project, institute=institute, filename=file, dataunits=dataunits, latunits=latunits, lonunits=lonunits, variable=variable, typeofvar=variable, typeofcal=caltype, varattribs=varattrib, globalattribs=attribs)


end

"""
    buildtimevec(str::String)Y

Construct the time vector from the netCDF file str

"""

function buildtimevec(str::String)

  # Time units
  ncI = NetCDF.ncinfo(str); # seems to be necessary. Otherwise can an inconsistent error when trying to load attributes
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
    leapDaysPer = sumleapyear(initDate, timeRaw[1] - 1)
    leapDaysPer2 = sumleapyear(initDate, timeRaw[end])
    startDate = initDate + Base.Dates.Day(convert(Int64, round(timeRaw[1]))) + Base.Dates.Day(leapDaysPer)
    endDate = initDate + Base.Dates.Day(convert(Int64, round(timeRaw[end]))) + Base.Dates.Day(leapDaysPer2)

    dateTmp = Date(startDate):Date(endDate)

    # REMOVE leap year (i.e. climate models use a no leap calendar)
    idx = (Dates.month.(dateTmp) .== 2) .&  (Dates.day.(dateTmp) .== 29)
    if length(idx) !== 1
      dateTmp = dateTmp[.!idx]
    else
      dateTmp = Array{Date}(dateTmp)
    end

elseif calType == "gregorian" || calType == "standard" || calType == "proleptic_gregorian"
    timeRaw = floor.(NetCDF.ncread(str, "time"))
    leapDaysPer = sumleapyear(initDate, timeRaw[1])
    leapDaysPer2 = sumleapyear(initDate, timeRaw[end])
    startDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[1])))
    endDate = initDate + Base.Dates.Day(convert(Int64,round(timeRaw[end])))
    dateTmp = Date(startDate):Date(endDate)

elseif calType == "360_day"

    timeRaw = floor.(NetCDF.ncread(str, "time"))
    leapDaysPer = sumleapyear(initDate, timeRaw[1] - 1)
    leapDaysPer2 = sumleapyear(initDate, timeRaw[end])

    startDate = initDate + Base.Dates.Day(convert(Int64, round(timeRaw[1]))) + Base.Dates.Day(leapDaysPer)
    endDate = initDate + Base.Dates.Day(convert(Int64, round(timeRaw[end]))) + Base.Dates.Day(leapDaysPer2)

    dateTmp = Date(startDate):Date(endDate)


  end
  # output date vector
  dateTmp = convert(Array{Date, 1}, dateTmp)
  return dateTmp
end

# function buildtimevec2(str::String)
#
#   # Time vectors
#   timevector = NetCDF.ncread(str, "time_vectors")
#   if size(timevector, 1) < size(timevector, 2)
#       timevector = timevector'
#   end
#
#   return Base.Date.(timevector)
# end

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

function corr_timevec(timeV, timefreq)

    if timefreq == "mon"
        year, month = Dates.year.(timeV), Dates.month.(timeV)
        years = unique(year)
        months = unique(month)
        timeout = Array{Date, 1}(length(years)*length(months))
        z = 1
        for iyear in years
            for imonth in months
                timeout[z] = Date(iyear, imonth)
                z += 1

            end
        end


    elseif timefreq == "day"
        error("have to do it")

    else
        error("wrong timefreq, needs to add it")

    end

    return timeout

end

"""
    timeresolution(timevec::Array{N,1} where N)

This function return the time resolution of the vector timevec, as obtained by the following call form the NetCDF.jl package

    timevec = NetCDF.ncread("netcdf_file.nc", "time")
"""

function timeresolution(timevec::Array{N,1} where N)

    # timevec = (NetCDF.ncread(str, "time"))
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

    function pr_timefactor(rez::String)

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


"""
    function spatialsubset(C::ClimGrid, poly::Array{N, 2})

Returns the spatial subset of ClimGrid C. The spatial subset is defined by a polygon poly

"""

function spatialsubset(C::ClimGrid, poly::Array{N, 2} where N)

    # Some checks for polygon poly

    if size(poly, 1) != 2 && size(poly, 2) == 2
        # transpose
        poly = poly'
    end


    lonsymbol = Symbol(C.dimension_dict["lon"])
    latsymbol = Symbol(C.dimension_dict["lat"])

    longrid = C.longrid
    latgrid = C.latgrid

    lon = C[1][Axis{lonsymbol}][:]
    lat = C[1][Axis{latsymbol}][:]

    msk = inpolygrid(longrid, latgrid, poly)
    idlon, idlat = findn(.!isnan.(msk))
    minXgrid = minimum(idlon)
    maxXgrid = maximum(idlon)
    minYgrid = minimum(idlat)
    maxYgrid = maximum(idlat)

    # Get DATA
    data = C[1].data

    if ndims(data) == 3
        data = data[:, minXgrid:maxXgrid, minYgrid:maxYgrid]
    elseif ndims(data) == 4
        data = data[:, minXgrid:maxXgrid, minYgrid:maxYgrid, :]
    end

    #new mask (e.g. representing the region of the polygon)
    msk = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]
    data = applymask(data, msk)

    # Get lon-lat for such region
    longrid = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
    latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
    lon = lon[minXgrid:maxXgrid]
    lat = lat[minYgrid:maxYgrid]

    dataOut = AxisArray(data, Axis{:time}(C[1][Axis{:time}][:]), Axis{lonsymbol}(lon), Axis{latsymbol}(lat))

    return ClimGrid(dataOut, longrid=longrid, latgrid=latgrid, msk=msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)

end

"""
    function temporalsubset(C::ClimGrid, start::Date, end::Date)

Returns the temporal subset of ClimGrid C. The temporal subset is defined by a start and end date.

"""

function temporalsubset(C::ClimGrid, startdate::Date, enddate::Date)

    # some checkups
    @argcheck startdate <= enddate
    @argcheck startdate >= C[1][Axis{:time}][1]
    @argcheck enddate <= C[1][Axis{:time}][end]

    # Temporal subset
    dataOut = C[1][Axis{:time}(startdate .. enddate)]

    return ClimGrid(dataOut, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)

end

function model_id(attrib::NCDatasets.Attributes)

    model = "N/A"

    try
        model = attrib["model_id"]
    catch
        try
            model = attrib["model"]
        end
    end

    return model
end

function experiment_id(attrib::NCDatasets.Attributes)

    experiment = "N/A"

    try
        experiment = attrib["experiment_id"]
    catch
        try
            experiment = attrib["experiment"]
        end
    end

    return experiment
end

function project_id(attrib::NCDatasets.Attributes)

    project = "N/A"

    try
        project = attrib["project_id"]
    catch
        try
            project = attrib["project"]
        end
    end

    return project
end

function institute_id(attrib::NCDatasets.Attributes)

    institute = "N/A"

    try
        institute = attrib["institute_id"]
    catch
        try
            institute = attrib["institute_id"]
        end
    end

    return institute

end

function frequency_var(attrib::NCDatasets.Attributes)

    frequency = "N/A"

    try
        frequency = attrib["frequency"]

    end

    return frequency
end

function runsim_id(attrib::NCDatasets.Attributes)

    runsim = "N/A"

    try
        runsim = attrib["parent_experiment_rip"]
    catch
        try
            runsim = attrib["driving_model_ensemble_member"]
        end
    end
    return runsim
end

function getdim_lat(ds::NCDatasets.Dataset)

    if sum(keys(ds.dim) .== "rlat") == 1
        return "rlat"
    elseif sum(keys(ds.dim) .== "lat") == 1
        return "lat"
    elseif sum(keys(ds.dim) .== "y") == 1
        return "y"
    else
        error("Manually verify x/lat dimension name")
    end

end

function getdim_lon(ds::NCDatasets.Dataset)

    if sum(keys(ds.dim) .== "rlon") == 1
        return "rlon"
    elseif sum(keys(ds.dim) .== "lon") == 1
        return "lon"
    elseif sum(keys(ds.dim) .== "x") == 1
        return "x"
    else
        error("Manually verify x/lat dimension name")
    end

end
