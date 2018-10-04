"""
    load(file::String, vari::String; poly = Array{Float64}([]), start_date::Tuple, end_date::Tuple, data_units::String = "")

Returns a ClimGrid type with the data in **file** of variable **vari** inside the polygon **poly**. Metadata is built-in the ClimGrid type, from the netCDF attributes.

Inside the ClimgGrid type, the data is stored into an AxisArray data type, with time, longitude/x and latitude/y dimensions.

The polygon provided should be in the -180, +180 longitude format. If the polygon crosses the International Date Line, the polygon should be splitted in multiple parts (i.e. multi-polygons).

Options for data_units are for precipitation : "mm", which converts the usual "kg m-2 s-1" unit found in netCDF files. For temperature : "Celsius", which converts the usual "Kelvin" unit.

Temporal subsetting can be done by providing start_date and end-date Tuples of length 1 (year), length 3 (year, month, day) or 6 (hour, minute, second).

**Note:** load uses [CF conventions](http://cfconventions.org/). If you are unable to read the netCDF file with load, the user will need to read it with low-level functions available in [NetCDF.jl package](https://github.com/JuliaGeo/NetCDF.jl) or [NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl) or re-create standartized netCDF files.
"""
function load(file::String, vari::String; poly = ([]), start_date::Tuple=(Inf,), end_date::Tuple=(Inf,), data_units::String = "")

  # TODO this file is a complete mess, but it works. Clean it up!
  # TODO create another function for orography extraction (2D dataset)

  # Get attributes
  ncI = NetCDF.ncinfo(file);
  attribs = NetCDF.ncinfo(file).gatts;
  ds = NCDatasets.Dataset(file)
  # ncI = NetCDF.ncinfo(file);
  attribs_dataset = ds.attrib

  # The following attributes should be set for netCDF files that follows CF conventions.
  # project, institute, model, experiment, frequency

  project = ClimateTools.project_id(attribs_dataset)
  institute = ClimateTools.institute_id(attribs_dataset)
  model = ClimateTools.model_id(attribs_dataset)
  experiment = ClimateTools.experiment_id(attribs_dataset)
  frequency = ClimateTools.frequency_var(attribs_dataset)
  runsim = ClimateTools.runsim_id(attribs_dataset)

  dataunits = ds[vari].attrib["units"]
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
  varattrib = Dict(ds[vari].attrib)

  if latname != "lat" && lonname != "lon" # means we don't have a "regular" grid
      latgrid = NetCDF.ncread(file, "lat")
      longrid = NetCDF.ncread(file, "lon")
      if @isdefined varattrib["grid_mapping"]
          map_dim = varattrib["grid_mapping"]
          map_attrib = Dict(ds[map_dim].attrib)
          map_attrib["grid_mapping"] = map_dim
      else
          error("File an issue on https://github.com/Balinus/ClimateTools.jl/issues to get the grid supported")
      end
  else # if no grid provided, create one
      longrid, latgrid = ndgrid(lon_raw, lat_raw)
      map_attrib = Dict(["grid_mapping_name" => "Regular_longitude_latitude"])
  end

  # =====================
  # TIME

  # Get time resolution

  rez = ClimateTools.timeresolution(NetCDF.ncread(file, "time"))

  # Construct time vector from info in netCDF file *str*
  timeV = buildtimevec(file, rez)
  if frequency == "mon"
      timeV = corr_timevec(timeV, frequency)
  end

  idxtimebeg, idxtimeend = timeindex(timeV, start_date, end_date, frequency)

  timeV = timeV[idxtimebeg:idxtimeend]

  # ==================
  # Spatial shift if grid is 0-360.
  rotatedgrid = false
  if sum(longrid .> 180) >= 1
      rotatedgrid = true

      # Shift 360 degrees grid to -180, +180 degrees
      longrid_flip = ClimateTools.shiftgrid_180_west_east(longrid)

      # Shift 360 degrees vector to -180, +180 degrees
      lon_raw_flip = ClimateTools.shiftvector_180_west_east(lon_raw)

  else
      longrid_flip = longrid # grid is already "flipped" by design
  end

  # ===================
  # GET DATA
  # data = ds[variable]
  data = NetCDF.open(file, vari)
  if !isempty(poly)

    # Test to see if the polygon crosses the meridian
    meridian = ClimateTools.meridian_check(poly)

    # Build mask based on provided polygon
    msk = inpolygrid(longrid_flip, latgrid, poly)

    if sum(isnan.(msk)) == length(msk) # no grid point insode polygon
        throw(error("No grid points found inside the provided polygon"))
    end

    if rotatedgrid
        # Regrid msk to original grid if the grid have been rotated to get proper index for data extraction
        msk = ClimateTools.shiftarray_east_west(msk, longrid_flip)
    end

    #Extract data based on mask
    data = ClimateTools.extractdata(data, msk, idxtimebeg, idxtimeend)

    #new mask (e.g. representing the region of the polygon)
    idlon, idlat = findn(.!isnan.(msk))
    minXgrid = minimum(idlon)
    maxXgrid = maximum(idlon)
    minYgrid = minimum(idlat)
    maxYgrid = maximum(idlat)
    msk = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]
    data = applymask(data, msk) # needed when polygon is not rectangular

    if rotatedgrid
        # Regrid msk to shifted grid if the grid have been rotated to get final mask
        #shiftarray_west_east(msk, longrid)
    end

    # Get lon_raw and lat_raw for such region
    lon_raw = lon_raw[minXgrid:maxXgrid]
    lat_raw = lat_raw[minYgrid:maxYgrid]

    if map_attrib["grid_mapping_name"] == "Regular_longitude_latitude"

        lon_raw = ClimateTools.shiftvector_180_west_east(lon_raw)

    end

    # Idem for longrid and latgrid
    if meridian
        longrid = ClimateTools.shiftgrid_180_east_west(longrid)#grideast, gridwest)
        longrid = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
        latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]

        # Re-translate to -180, 180 if 0, 360

        if rotatedgrid

            longrid_flip = ClimateTools.shiftgrid_180_west_east(longrid)

            data = permute_west_east(data, longrid)#idxwest, idxeast)
            msk = ClimateTools.permute_west_east(msk, longrid)

            # TODO Try to trim padding when meridian is crossed and model was on a 0-360 coords
            # idlon, idlat = findn(.!isnan.(msk))
            # minXgrid = minimum(idlon)
            # maxXgrid = maximum(idlon)
            # minYgrid = minimum(idlat)
            # maxYgrid = maximum(idlat)
            #
            # msk = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]
            # data = data[:, minXgrid:maxXgrid, minYgrid:maxYgrid]
            # lon_raw_flip = lon_raw[minXgrid:maxXgrid]
            #
            # longrid_flip = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
            # longrid_flip = ClimateTools.shiftgrid_180_west_east(longrid_flip)
            # latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]

            #     data = applymask(data, msk)



        end




    else

        if rotatedgrid # flip in original grid
            longrid = shiftgrid_180_east_west(longrid) #grideast, gridwest)
        end

        longrid = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
        latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]


        if rotatedgrid

            longrid_flip = shiftgrid_180_west_east(longrid)

            lon_raw_flip = shiftvector_180_east_west(lon_raw)

            data = permute_west_east(data, longrid)#idxwest, idxeast)
            msk = ClimateTools.permute_west_east(msk, longrid)

        end

    end

  elseif isempty(poly) # no polygon clipping
      msk = Array{Float64}(ones((size(data, 1), size(data, 2))))
    # if ndims(data) == 3
      data = extractdata(data, msk, idxtimebeg, idxtimeend)
    # elseif ndims(data) == 4
        # data = extractdata(data, msk, idxtimebeg, idxtimeend)
    # end


    if rotatedgrid
        # Flip data "west-east"
        data = permute_west_east(data, longrid)
    end

  end

  # # Replace fillvalues with NaN
  fillval = NetCDF.ncgetatt(file, vari, "_FillValue")
  data[data .== fillval] = NaN

  if rotatedgrid
      longrid .= longrid_flip
      lon_raw .= lon_raw_flip
  end

  # Convert units of optional argument data_units is provided
  if data_units == "Celsius" && (vari == "tas" || vari == "tasmax" || vari == "tasmin") && dataunits == "K"
    data = data .- 273.15
    dataunits = "Celsius"
  end

  if data_units == "mm" && vari == "pr" && (dataunits == "kg m-2 s-1" || dataunits == "mm s-1")

    rez = timeresolution(NetCDF.ncread(file, "time"))
    factor = pr_timefactor(rez)
    data = data .* factor
    if rez != "N/A"
        dataunits = string("mm/",rez)
    else
        dataunits = "mm"
    end
    varattrib["standard_name"] = "precipitation"
  end

  # Create AxisArray from variable "data"
  if ndims(data) == 3
    # Convert data to AxisArray
    dataOut = AxisArray(data, Axis{Symbol(lonname)}(lon_raw), Axis{Symbol(latname)}(lat_raw), Axis{:time}(timeV))
  elseif ndims(data) == 4 # this imply a 3D field (height component)
    # Get level vector
    plev = ds["plev"][:]#NetCDF.ncread(file, "plev")
    # Convert data to AxisArray
    dataOut = AxisArray(data, Axis{Symbol(lonname)}(lon_raw), Axis{Symbol(latname)}(lat_raw), Axis{:plev}(plev), Axis{:time}(timeV))
  else
    throw(error("load takes only 3D and 4D variables for the moment"))
  end

  close(ds)
  NetCDF.ncclose(file)

  return ClimGrid(dataOut, longrid=longrid, latgrid=latgrid, msk=msk, grid_mapping=map_attrib, dimension_dict=dimension_dict, model=model, frequency=frequency, experiment=experiment, run=runsim, project=project, institute=institute, filename=file, dataunits=dataunits, latunits=latunits, lonunits=lonunits, variable=vari, typeofvar=vari, typeofcal=caltype, varattribs=varattrib, globalattribs=attribs)


end

"""
    load(files::Array{String,1}, vari::String; poly = ([]), start_date::Date = Date(-4000), end_date::Date = Date(-4000), data_units::String = "")

Loads and merge the files contained in the arrar files.
"""
function load(files::Array{String,1}, vari::String; poly = ([]), start_date::Tuple=(Inf,), end_date::Tuple=(Inf,), data_units::String = "")

    C = [] # initialize # TODO better initialization
    nfiles = length(files)

    p = Progress(nfiles, 3, "Loading files: ")

    for ifile = 1:nfiles

        datatmp = load(files[ifile], vari, poly = poly, start_date=start_date, end_date=end_date, data_units=data_units)

        if ifile == 1
            C = datatmp
        else
            C = merge(C, datatmp)
        end
        next!(p)
    end
    return C
end


"""
    load2D(file::String, vari::String; poly=[], data_units::String="")

Returns a 2D array. Should be used for *fixed* data, such as orography
"""
function load2D(file::String, vari::String; poly=[], data_units::String="")
    # Get attributes
    ncI = NetCDF.ncinfo(file);
    attribs = NetCDF.ncinfo(file).gatts;
    ds = NCDatasets.Dataset(file)
    # ncI = NetCDF.ncinfo(file);
    attribs_dataset = ds.attrib

    project = ClimateTools.project_id(attribs_dataset)
    institute = ClimateTools.institute_id(attribs_dataset)
    model = ClimateTools.model_id(attribs_dataset)
    experiment = ClimateTools.experiment_id(attribs_dataset)
    frequency = ClimateTools.frequency_var(attribs_dataset)
    runsim = ClimateTools.runsim_id(attribs_dataset)

    dataunits = ds[vari].attrib["units"]
    latunits = ds["lat"].attrib["units"]
    lonunits = ds["lon"].attrib["units"]

    # Get dimensions names
    latname = getdim_lat(ds)
    lonname = getdim_lon(ds)

    # Create dict with latname and lonname
    dimension_dict = Dict(["lon" => lonname, "lat" => latname])

    lat_raw = NetCDF.ncread(file, latname)
    lon_raw = NetCDF.ncread(file, lonname)

    # Get variable attributes
    varattrib = Dict(ds[vari].attrib)

    if latname != "lat" && lonname != "lon" # means we don't have a "regular" grid
        latgrid = NetCDF.ncread(file, "lat")
        longrid = NetCDF.ncread(file, "lon")
        if ClimateTools.@isdefined varattrib["grid_mapping"]
            map_dim = varattrib["grid_mapping"]
            map_attrib = Dict(ds[map_dim].attrib)
            map_attrib["grid_mapping"] = map_dim
        else
            error("File an issue on https://github.com/Balinus/ClimateTools.jl/issues to get the grid supported")
        end
    else # if no grid provided, create one
        longrid, latgrid = ndgrid(lon_raw, lat_raw)
        map_attrib = Dict(["grid_mapping_name" => "Regular_longitude_latitude"])
    end

    # =====================
    # TIME

    # ==================
    # Spatial shift if grid is 0-360.
    rotatedgrid = false
    if sum(longrid .> 180) >= 1
        rotatedgrid = true

        # Shift 360 degrees grid to -180, +180 degrees
        longrid_flip = ClimateTools.shiftgrid_180_west_east(longrid)

        # Shift 360 degrees vector to -180, +180 degrees
        lon_raw_flip = ClimateTools.shiftvector_180_west_east(lon_raw)

    else
        longrid_flip = longrid # grid is already "flipped" by design
    end

    # ===================
    # GET DATA
    # data = ds[variable]
    data = NetCDF.open(file, vari)
    if !isempty(poly)

      # Test to see if the polygon crosses the meridian
      meridian = ClimateTools.meridian_check(poly)

      # Build mask based on provided polygon
      msk = inpolygrid(longrid_flip, latgrid, poly)

      if sum(isnan.(msk)) == length(msk) # no grid point insode polygon
          throw(error("No grid points found inside the provided polygon"))
      end

      if rotatedgrid
          # Regrid msk to original grid if the grid have been rotated to get proper index for data extraction
          msk = ClimateTools.shiftarray_east_west(msk, longrid_flip)
      end

      #Extract data based on mask
      data = ClimateTools.extractdata2D(data, msk)

      #new mask (e.g. representing the region of the polygon)
      # idlon, idlat = findn(.!isnan.(msk)) # DEPRECATED SEE NEXT "begin...end"
      
      begin
        I = Base.findall(!isnan, msk)
        idlon, idlat = (getindex.(I, 1), getindex.(I, 2))
      end

      minXgrid = minimum(idlon)
      maxXgrid = maximum(idlon)
      minYgrid = minimum(idlat)
      maxYgrid = maximum(idlat)
      msk = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]
      data = applymask(data, msk) # needed when polygon is not rectangular

      if rotatedgrid
          # Regrid msk to shifted grid if the grid have been rotated to get final mask
          #shiftarray_west_east(msk, longrid)
      end

      # Get lon_raw and lat_raw for such region
      lon_raw = lon_raw[minXgrid:maxXgrid]
      lat_raw = lat_raw[minYgrid:maxYgrid]

      if map_attrib["grid_mapping_name"] == "Regular_longitude_latitude"

          lon_raw = ClimateTools.shiftvector_180_west_east(lon_raw)

      end

      # Idem for longrid and latgrid
      if meridian
          longrid = ClimateTools.shiftgrid_180_east_west(longrid)#grideast, gridwest)
          longrid = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
          latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]

          # Re-translate to -180, 180 if 0, 360

          if rotatedgrid

              longrid_flip = ClimateTools.shiftgrid_180_west_east(longrid)

              data = permute_west_east(data, longrid)#idxwest, idxeast)
              msk = ClimateTools.permute_west_east(msk, longrid)

              # TODO Try to trim padding when meridian is crossed and model was on a 0-360 coords
              # idlon, idlat = findn(.!isnan.(msk))
              # minXgrid = minimum(idlon)
              # maxXgrid = maximum(idlon)
              # minYgrid = minimum(idlat)
              # maxYgrid = maximum(idlat)
              #
              # msk = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]
              # data = data[:, minXgrid:maxXgrid, minYgrid:maxYgrid]
              # lon_raw_flip = lon_raw[minXgrid:maxXgrid]
              #
              # longrid_flip = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
              # longrid_flip = ClimateTools.shiftgrid_180_west_east(longrid_flip)
              # latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]

              #     data = applymask(data, msk)



          end




      else

          if rotatedgrid # flip in original grid
              longrid = shiftgrid_180_east_west(longrid) #grideast, gridwest)
          end

          longrid = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
          latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]


          if rotatedgrid

              longrid_flip = shiftgrid_180_west_east(longrid)

              lon_raw_flip = shiftvector_180_east_west(lon_raw)

              data = permute_west_east(data, longrid)#idxwest, idxeast)
              msk = ClimateTools.permute_west_east(msk, longrid)

          end

      end

    elseif isempty(poly) # no polygon clipping
        msk = Array{Float64}(ones((size(data, 1), size(data, 2))))
      # if ndims(data) == 3
        data = extractdata2D(data, msk)
      # elseif ndims(data) == 4
          # data = extractdata(data, msk, idxtimebeg, idxtimeend)
      # end


      if rotatedgrid
          # Flip data "west-east"
          data = permute_west_east(data, longrid)
      end

    end

    # # Replace fillvalues with NaN
    fillval = NetCDF.ncgetatt(file, vari, "_FillValue")
    data[data .== fillval] = NaN

    if rotatedgrid
        longrid = longrid_flip
        lon_raw = lon_raw_flip
    end

    # Convert units of optional argument data_units is provided
    if data_units == "Celsius" && (vari == "tas" || vari == "tasmax" || vari == "tasmin") && dataunits == "K"
      data = data - 273.15
      dataunits = "Celsius"
    end

    if data_units == "mm" && vari == "pr" && (dataunits == "kg m-2 s-1" || dataunits == "mm s-1")

      rez = timeresolution(NetCDF.ncread(file, "time"))
      factor = pr_timefactor(rez)
      data = data .* factor
      if rez != "N/A"
          dataunits = string("mm/",rez)
      else
          dataunits = "mm"
      end
      varattrib["standard_name"] = "precipitation"
    end

    # Create AxisArray from variable "data"

    # Convert data to AxisArray
    dataOut = AxisArray(data, Axis{Symbol(lonname)}(lon_raw), Axis{Symbol(latname)}(lat_raw))


    close(ds)
    NetCDF.ncclose(file)

    return ClimGrid(dataOut, longrid=longrid, latgrid=latgrid, msk=msk, grid_mapping=map_attrib, dimension_dict=dimension_dict, model=model, frequency=frequency, experiment=experiment, run=runsim, project=project, institute=institute, filename=file, dataunits=dataunits, latunits=latunits, lonunits=lonunits, variable=vari, typeofvar=vari, typeofcal="fixed", varattribs=varattrib, globalattribs=attribs)




end


"""
    buildtimevec(str::String)

Construct the time vector from the netCDF file str

"""
function buildtimevec(str::String, rez)

  # Time units
  ncI = NetCDF.ncinfo(str); # seems to be necessary. Otherwise can an inconsistent error when trying to load attributes
  units = NetCDF.ncgetatt(str, "time", "units") # get starting date
  # m = match(r"(\d+)[-.\/](\d+)[-.\/](\d+)", units, 1)
  m = match(r"(\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})", units, 1) # match a date from string
  if isempty(fieldnames(typeof(m)))#@isdefined m.captures
      m = match(r"(\d+)[-.\/](\d+)[-.\/](\d+)", units, 1)
  end

  arrdate = parse.(Int, m.captures)
  if length(arrdate) == 6
      initDate = DateTime(arrdate[1], arrdate[2], arrdate[3], arrdate[4], arrdate[5], arrdate[6])
  elseif length(arrdate) == 3
      initDate = DateTime(arrdate[1], arrdate[2], arrdate[3])
  end
  # daysfrom = m.match # get only the date ()"yyyy-mm-dd" format)
  # initDate = Date(daysfrom, "yyyy-mm-dd")

  period = getperiod(rez)
  timeRaw = NetCDF.ncread(str, "time")

  # Calendar type
  calType = NetCDF.ncgetatt(str, "time", "calendar")
  if calType == "noleap" || calType == "365_day"
    nDays = 365
    # get time of netCDF file *str*
    # timeRaw = floor.(NetCDF.ncread(str, "time"))

    leapDaysPer = sumleapyear(initDate, timeRaw[1] - 1)
    leapDaysPer2 = sumleapyear(initDate, timeRaw[end])
    startDate = initDate + Dates.Day(convert(Int64, floor(timeRaw[1]))) + Dates.Day(leapDaysPer)
    endDate = initDate + Dates.Day(convert(Int64, ceil(timeRaw[end]))) + Dates.Day(leapDaysPer2) - period

    # period = getperiod(rez)

    dateTmp = DateTime(startDate):period:DateTime(endDate)

    # REMOVE leap year (i.e. climate models use a no leap calendar)
    idx = (Dates.month.(dateTmp) .== 2) .&  (Dates.day.(dateTmp) .== 29)
    if length(idx) !== 1
      dateTmp = dateTmp[.!idx]
    else
      dateTmp = Array{Date}(dateTmp)
    end

elseif calType == "gregorian" || calType == "standard" || calType == "proleptic_gregorian"
    # timeRaw = floor.(NetCDF.ncread(str, "time"))
    leapDaysPer = sumleapyear(initDate, timeRaw[1])
    leapDaysPer2 = sumleapyear(initDate, timeRaw[end])

    if typeof(timeRaw[1]) == Int8 || typeof(timeRaw[1]) == Int16 || typeof(timeRaw[1]) == Int32 || typeof(timeRaw[1]) == Int64

        startDate = initDate + Dates.Day(floor(timeRaw[1]))
        endDate = initDate + Dates.Day(floor(timeRaw[end]))
    else
        startDate = initDate + Dates.Day(convert(Int64,floor(timeRaw[1])))
        endDate = initDate + Dates.Day(convert(Int64,floor(timeRaw[end])))# - period
    end
    dateTmp = DateTime(startDate):period:DateTime(endDate)

elseif calType == "360_day"
    throw(error("360_day type of calendar not yet supported"))

    # timeRaw = floor.(NetCDF.ncread(str, "time"))
    leapDaysPer = sumleapyear(initDate, timeRaw[1] - 1)
    leapDaysPer2 = sumleapyear(initDate, timeRaw[end])

    startDate = initDate + Dates.Day(convert(Int64, round(timeRaw[1]))) + Dates.Day(leapDaysPer)
    endDate = initDate + Dates.Day(convert(Int64, round(timeRaw[end]))) + Dates.Day(leapDaysPer2)

    dateTmp = DateTime(startDate):period:DateTime(endDate)


  end
  # output date vector
  dateTmp = convert(Array{DateTime, 1}, dateTmp)
  return dateTmp
end

"""
    getperiod(rez::String)

Returns the Dates resolution (e.g. Dates.Hour(1) for hourly data, Dates.Hour(3) for 3-hours data, etc..).
"""
function getperiod(rez::String)

    if rez == "1h"
        period = Dates.Hour(1)
    elseif rez == "3h"
        period = Dates.Hour(3)
    elseif rez == "6h"
        period = Dates.Hour(6)
    elseif rez == "12h"
        period = Dates.Hour(12)
    elseif rez == "24h"
        period = Dates.Hour(24)
    elseif rez == "Yearly"
        period = Dates.Year(1)
    else
        period = Dates.Hour(24)
    end

    return period
end


"""
Number of leap years in date vector

    sumleapyear(dates::StepRange{Date, Dates.Day})

    sumleapyear(initDate::Date, timeRaw)
"""
function sumleapyear(initDate::Dates.TimeType, timeRaw)

  out = 0::Int
  endDate = initDate + Dates.Day(convert(Int64,round(timeRaw[1])))
  years = unique(Dates.year.(initDate:Day(1):endDate))
  # Sum over time vector
  for idx = 1:length(years)
    if Dates.isleapyear(years[idx])
      out += 1
    end

  end

  return out

end

function sumleapyear(dates::StepRange{Date, Dates.Day})

  out = 0::Int
  endDate = dates[end]#initDate + Dates.Day(convert(Int64,round(timeRaw[1])))
  years = unique(Dates.year.(dates))
  # Sum over time vector
  for idx = 1:length(years)
    if Dates.isleapyear(years[idx])
      out += 1
    end
  end

  return out

end

"""
    shapefile_coords(poly::Shapefile.Polygon)

This function return the polygons contained in shp.shapes[i]. It returns the x and y coordinates vectors.

See also [`shapefile_coords_poly`](@ref), which returns a polygon that ca be used for data extraction of the [`load`](@ref).

"""
function shapefile_coords(poly::Shapefile.Polygon)
    start_indices = poly.parts .+ 1
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
    shapefile_coords_poly(poly::Shapefile.Polygon)

Return the polygons contained in shp.shapes[i]. It returns an array containing the polygons.

See also [`shapefile_coords`](@ref), which returns vectors as opposed to array. Returned polygon is consistent with the data extraction of the [`load`](@ref) function.

"""
function shapefile_coords_poly(poly::Shapefile.Polygon)
    x, y = shapefile_coords(poly)
    return [x y]'
end

"""
    extractpoly(file::String, n::Int)

Returns the n-th polygon contained in *file*.
"""
function extractpoly(file::String, n::Int)

    shp = open(file) do fd
        read(fd, Shapefile.Handle)
    end

    poly = shapefile_coords_poly(shp.shapes[n])
    return poly
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

Return the time resolution of the vector timevec.

"""
function timeresolution(timevec::Array{N,1} where N)

    # timevec = (NetCDF.ncread(str, "time"))
    if length(timevec) > 1
        timediff = diff(timevec)[2]
        if timediff == 1. || timediff == 1
            return "24h"
        elseif round(timediff, digits=5) == round(12/24, digits=5)
            return "12h"
        elseif round(timediff, digits=5) == round(6/24, digits=5)
            return "6h"
        elseif round(timediff, digits=5) == round(3/24, digits=5)
            return "3h"
        elseif round(timediff, digits=5) == round(1/24, digits=5)
            return "1h"
        elseif round(timediff, digits=5) == 365.0 || round(timediff, digits=5) == 366.0 || round(timediff, digits=5) == 360.0
            return "Yearly"
        end
    else
        return "N/A"
    end
end

"""
    function pr_timefactor(rez::String)

Return the time factor that should be applied to precipitation to get accumulation for resolution "rez"

"""
function pr_timefactor(rez::String)

    if rez == "24h"
        return 86400.0
    elseif rez == "12h"
        return 43200.0
    elseif rez == "6h"
        return 21600.0
    elseif rez == "3h"
        return 10800.0
    elseif rez == "1h"
        return 3600.0
    elseif rez == "N/A"
        return 1.0
    end

end


"""
    spatialsubset(C::ClimGrid, poly::Array{N, 2})

Returns the spatial subset of ClimGrid C. The spatial subset is defined by the polygon poly, defined on a -180, +180 longitude reference.

"""
function spatialsubset(C::ClimGrid, poly)

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
    # idlon, idlat = findn(.!isnan.(msk))
    begin
        I = findall(!isnan, msk)
        idlon, idlat = (getindex.(I, 1), getindex.(I, 2))
    end
    minXgrid = minimum(idlon)
    maxXgrid = maximum(idlon)
    minYgrid = minimum(idlat)
    maxYgrid = maximum(idlat)

    # Get DATA
    data = C[1].data

    if ndims(data) == 3
        data = data[minXgrid:maxXgrid, minYgrid:maxYgrid, :]
    elseif ndims(data) == 4
        data = data[minXgrid:maxXgrid, minYgrid:maxYgrid, :, :]
    end

    #new mask (e.g. representing the region of the polygon)
    msk = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]
    data = applymask(data, msk)

    # Get lon-lat for such region
    longrid = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
    latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
    lon = lon[minXgrid:maxXgrid]
    lat = lat[minYgrid:maxYgrid]

    if ndims(data) == 3
      dataOut = AxisArray(data, Axis{lonsymbol}(lon), Axis{latsymbol}(lat),  Axis{:time}(C[1][Axis{:time}][:]))
    elseif ndims(data) == 4
      dataOut = AxisArray(data, Axis{lonsymbol}(lon), Axis{latsymbol}(lat), Axis{:level}(C[1][Axis{:plev}][:]), Axis{:time}(C[1][Axis{:time}][:]))
    end

    return ClimGrid(dataOut, longrid=longrid, latgrid=latgrid, msk=msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)

end

"""
    function temporalsubset(C::ClimGrid, startdate::Date, enddate::Date)

Returns the temporal subset of ClimGrid C. The temporal subset is defined by a start and end date.

"""
function temporalsubset(C::ClimGrid, datebeg::Tuple, dateend::Tuple)

    startdate = buildtimetype(datebeg)
    enddate = buildtimetype(dateend)

    # some checkups
    @argcheck startdate <= enddate
    @argcheck startdate >= C[1][Axis{:time}][1]
    @argcheck enddate <= C[1][Axis{:time}][end]

    # Temporal subset
    dataOut = C[1][Axis{:time}(startdate .. enddate)]

    return ClimGrid(dataOut, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)

end

"""
    periodsubset(C::ClimGrid, startmonth::Int64, endmonth::Int64)

Return the temporal subset of ClimGrid C based on months.
"""
function periodsubset(C::ClimGrid, startmonth::Int64, endmonth::Int64)

    @argcheck startmonth >= minimum(Dates.month.(C[1][Axis{:time}][:]))
    @argcheck startmonth <= maximum(Dates.month.(C[1][Axis{:time}][:]))

    if startmonth <= endmonth
        # Each matrix [:,:,i] represent data for a day
        datain = C.data.data
        # Date vector
        datevecin = C[1][Axis{:time}][:]
        # Where are the data between startmonth and endmonth
        index = (Dates.month.(datevecin) .<= endmonth) .&  (Dates.month.(datevecin) .>= startmonth)
        # Keep only data between startmonth and endmonth
        dataout = datain[:,:,index]
        datevecout = datevecin[index]
        # Create the ClimGrid output
        lonsymbol = Symbol(C.dimension_dict["lon"])
        latsymbol = Symbol(C.dimension_dict["lat"])
        axisout = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}][:]), Axis{latsymbol}(C[1][Axis{latsymbol}][:]),    Axis{:time}(datevecout))

    elseif endmonth < startmonth
        # Each matrix [:,:,i] represent data for a day
        datain = C.data.data
        # Date vector
        datevecin = C[1][Axis{:time}][:]
        # Where are the data between startmonth and endmonth
        index = (Dates.month.(datevecin) .<= endmonth) .|  (Dates.month.(datevecin) .>= startmonth)
        # Keep only data between startmonth and endmonth
        dataout = datain[:,:,index]
        datevecout = datevecin[index]
        # Create the ClimGrid output
        lonsymbol = Symbol(C.dimension_dict["lon"])
        latsymbol = Symbol(C.dimension_dict["lat"])
        axisout = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}][:]), Axis{latsymbol}(C[1][Axis{latsymbol}][:]), Axis{:time}(datevecout))
    end

    return ClimGrid(axisout, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping,dimension_dict=C.dimension_dict, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project,institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable,typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

"""
    periodsubset(C::ClimGrid, season::String)

Return the temporal subset of ClimGrid C for a given season. Season options are: "DJF" (December-February), "MAM" (March-May), "JJA" (June-August), "SON" (September-November)
"""
function periodsubset(C::ClimGrid, season::String)
    if lowercase(season) == "djf"
      D = periodsubset(C, 12, 2)
    elseif lowercase(season) == "mam"
      D = periodsubset(C, 3, 5)
    elseif lowercase(season) == "jja"
      D = periodsubset(C, 6, 8)
    elseif lowercase(season) == "son"
      D = periodsubset(C, 9, 11)
    else
      error("Wrong season name. Options are DJF (December-February; Winter), MAM (March-May; Spring), JJA (June-August; Summer) and SON (September-November; Fall)")
    end

    return D
end

model_id(attrib::NCDatasets.Attributes) = get(attrib,"model_id",get(attrib,"model","N/A"))

experiment_id(attrib::NCDatasets.Attributes) = get(attrib,"experiment_id",get(attrib,"experiment","N/A"))

project_id(attrib::NCDatasets.Attributes) = get(attrib,"project_id",get(attrib,"project","N/A"))

institute_id(attrib::NCDatasets.Attributes) = get(attrib,"institute_id",get(attrib,"institute","N/A"))

frequency_var(attrib::NCDatasets.Attributes) = get(attrib,"frequency","N/A")

runsim_id(attrib::NCDatasets.Attributes) = get(attrib, "parent_experiment_rip", get(attrib,"driving_model_ensemble_member","N/A"))


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

function shiftgrid_180_west_east(longrid)

    for ilon in eachindex(longrid)
        if longrid[ilon] >= 180
            longrid[ilon] = longrid[ilon] - 360
        end
    end

    ieast = longrid .>= 0
    iwest = longrid .< 0

    grideast = reshape(longrid[ieast], :, size(longrid, 2))
    gridwest = reshape(longrid[iwest], :, size(longrid, 2))

    longrid_flip = vcat(gridwest, grideast)

    return longrid_flip

end


function shiftgrid_180_east_west(longrid)

    for ilon in eachindex(longrid)
        if longrid[ilon] >= 180
            longrid[ilon] = longrid[ilon] - 360
        end
    end

    ieast = longrid .>= 0
    iwest = longrid .< 0

    grideast = reshape(longrid[ieast], :, size(longrid, 2))
    gridwest = reshape(longrid[iwest], :, size(longrid, 2))

    longrid_flip = vcat(grideast, gridwest)

    return longrid_flip

end

function shiftarray_west_east(data, longrid_flip)

    ieast = longrid_flip .>= 0
    iwest = longrid_flip .< 0

    msk = permute_west_east2D(data, iwest, ieast)
    return msk

end

function shiftarray_east_west(data, longrid_flip)

    ieast = longrid_flip .>= 0
    iwest = longrid_flip .< 0

    msk = permute_east_west2D(data, iwest, ieast)
    return msk

end


function shiftvector_180_east_west(lon_raw)

    for ilon in eachindex(lon_raw)
        if lon_raw[ilon] >= 180
            lon_raw[ilon] = lon_raw[ilon] - 360
        end
    end
    ieast_vec = lon_raw .>= 0
    iwest_vec = lon_raw .< 0
    lon_raw_flip = [lon_raw[ieast_vec];lon_raw[iwest_vec]]
    return lon_raw_flip

end

function shiftvector_180_west_east(lon_raw)

    for ilon in eachindex(lon_raw)
        if lon_raw[ilon] >= 180
            lon_raw[ilon] = lon_raw[ilon] - 360
        end
    end

    ieast_vec = lon_raw .>= 0
    iwest_vec = lon_raw .< 0

    lon_raw_flip = [lon_raw[iwest_vec];lon_raw[ieast_vec]]

    return lon_raw_flip

end

"""
    extractdata(data, msk, idxtimebeg, idxtimeend)

Returns the data contained in netCDF file, using the appropriate mask and time index. Used internally by `load`.

"""
function extractdata(data, msk, idxtimebeg, idxtimeend)

    # idlon, idlat = findn(.!isnan.(msk))
    begin
        I = Base.findall(!isnan, msk)
        idlon, idlat = (getindex.(I, 1), getindex.(I, 2))
    end
    minXgrid = minimum(idlon)
    maxXgrid = maximum(idlon)
    minYgrid = minimum(idlat)
    maxYgrid = maximum(idlat)

    if ndims(data) == 3
        data = data[minXgrid:maxXgrid, minYgrid:maxYgrid, idxtimebeg:idxtimeend]
        # Permute dims
        # data = permutedims(data, [3, 1, 2])
    elseif ndims(data) == 4
        data = data[minXgrid:maxXgrid, minYgrid:maxYgrid, :, idxtimebeg:idxtimeend]
        # Permute dims
        # data = permutedims(data, [4, 1, 2, 3])
    end

    return data
end

function extractdata2D(data, msk)

    # idlon, idlat = findn(.!isnan.(msk))
    begin
        I = Base.findall(!isnan, msk)
        idlon, idlat = (getindex.(I, 1), getindex.(I, 2))
    end
    minXgrid = minimum(idlon)
    maxXgrid = maximum(idlon)
    minYgrid = minimum(idlat)
    maxYgrid = maximum(idlat)

    data = data[minXgrid:maxXgrid, minYgrid:maxYgrid]

    return data

end

function meridian_check(poly)

    minpoly = minimum(poly[1, .!isnan.(poly[1, :])])
    maxpoly = maximum(poly[1, .!isnan.(poly[1, :])])

    meridian = false
    if sign(minpoly)*sign(maxpoly) == sign(-1.0) # polygon crosses the meridian
       meridian = true
    end

    return meridian

end
