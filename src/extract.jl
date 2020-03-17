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

  # Get attributes
  ds = NCDatasets.Dataset(file)
  attribs_dataset = ds.attrib
  attribs = Dict(attribs_dataset)

  # Data pointer
  data_pointer = ds[vari]

  # The following attributes should be set for netCDF files that follows CF conventions.
  # project, institute, model, experiment, frequency
  project = ClimateTools.project_id(attribs_dataset)
  institute = ClimateTools.institute_id(attribs_dataset)
  model = ClimateTools.model_id(attribs_dataset)
  experiment = ClimateTools.experiment_id(attribs_dataset)
  frequency = ClimateTools.frequency_var(attribs_dataset)
  runsim = ClimateTools.runsim_id(attribs_dataset)
  grid_mapping = ClimateTools.get_mapping(keys(ds))#, vari)

  if grid_mapping == "Regular_longitude_latitude"
      latstatus = false
  else
      latstatus = true
  end

  # Get dimensions names
  lonname = get_dimname(ds, "X")
  latname = get_dimname(ds, "Y")
  timename = get_dimname(ds, "T")
  if ndims(data_pointer) == 4
      levname = get_dimname(ds, "Z")
      levunits = ds[levname].attrib["units"]
  end
  # Dimension vector
  lat_raw = nomissing(ds[latname][:], NaN)
  lon_raw = nomissing(ds[lonname][:], NaN)

  # Create dict with latname and lonname
  if ndims(data_pointer) == 3
      dimension_dict = Dict(["lon" => lonname,
                             "lat" => latname,
                             "time" => timename])

  elseif ndims(data_pointer) == 4
      dimension_dict = Dict(["lon" => lonname,
                             "lat" => latname,
                             "height" => levname,
                             "time" => "time"])
  end

  # Get units
  dataunits = ds[vari].attrib["units"]
  latunits = ds[latname].attrib["units"]
  lonunits = ds[lonname].attrib["units"]

  # Calendar
  caltype = get_calendar(ds, timename)

  # Get variable attributes
  varattrib = Dict(ds[vari].attrib)

  if latstatus # means we don't have a "regular" grid
      # Get names of grid
      latgrid_name = ClimateTools.latgridname(ds)
      longrid_name = ClimateTools.longridname(ds)

      latgrid = nomissing(ds[latgrid_name][:], NaN)
      longrid = nomissing(ds[longrid_name][:], NaN)

      # Ensure we have a grid
      if ndims(latgrid) == 1 && ndims(longrid) == 1
          longrid, latgrid = ndgrid(lon_raw, lat_raw)
          map_attrib = Dict(["grid_mapping" => "Regular_longitude_latitude"])
      end

      map_attrib = build_grid_mapping(ds, grid_mapping)
      varattrib["grid_mapping"] = grid_mapping

  else # if no grid provided, create one
      longrid, latgrid = ndgrid(lon_raw, lat_raw)
      map_attrib = Dict(["grid_mapping" => grid_mapping])
      varattrib["grid_mapping"] = grid_mapping
  end

  # =====================
  # TIME
  # Get time resolution
  timeV = ds[timename][:]
  if frequency == "N/A" || !ClimateTools.@isdefined frequency
      try
          try
              frequency = string(diff(timeV)[2])
          catch
              frequency = string(diff(timeV)[1])
          end
      catch
          frequency = "N/A"
      end
  end

  timeattrib = Dict(ds[timename].attrib)
  T = typeof(timeV[1])
  idxtimebeg, idxtimeend = timeindex(timeV, start_date, end_date, T)
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
  if !isempty(poly)

    # Test to see if the polygon crosses the meridian
    meridian = ClimateTools.meridian_check(poly)

    # Build mask based on provided polygon
    msk = inpolygrid(longrid_flip, latgrid, poly)

    if sum(isnan.(msk)) == length(msk) # no grid point inside polygon
        throw(error("No grid points found inside the provided polygon"))
    end

    if rotatedgrid
        # Regrid msk to original grid if the grid have been rotated to get proper index for data extraction
        msk = ClimateTools.shiftarray_east_west(msk, longrid_flip)
    end

    #Extract data based on mask
    data_ext = ClimateTools.extractdata(data_pointer, msk, idxtimebeg, idxtimeend)
    data_ext = nomissing(data_ext, NaN)


    #new mask (e.g. representing the region of the polygon)
    begin
        I = Base.findall(!isnan, msk)
        idlon, idlat = (getindex.(I, 1), getindex.(I, 2))
      end
    minXgrid = minimum(idlon)
    maxXgrid = maximum(idlon)
    minYgrid = minimum(idlat)
    maxYgrid = maximum(idlat)
    msk = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]
    data_mask = applymask(data_ext, msk) # needed when polygon is not rectangular

    if rotatedgrid
        # Regrid msk to shifted grid if the grid have been rotated to get final mask
        #shiftarray_west_east(msk, longrid)
    end

    # Get lon_raw and lat_raw for such region
    lon_raw = lon_raw[minXgrid:maxXgrid]
    lat_raw = lat_raw[minYgrid:maxYgrid]

    if map_attrib["grid_mapping"] == "Regular_longitude_latitude"

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

            data = permute_west_east(data_mask, longrid)#idxwest, idxeast)
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

        else
            data = data_mask
        end

    else

        if rotatedgrid # flip in original grid
            longrid = ClimateTools.shiftgrid_180_east_west(longrid) #grideast, gridwest)
        end

        longrid = longrid[minXgrid:maxXgrid, minYgrid:maxYgrid]
        latgrid = latgrid[minXgrid:maxXgrid, minYgrid:maxYgrid]

        if rotatedgrid

            longrid_flip = shiftgrid_180_west_east(longrid)

            lon_raw_flip = shiftvector_180_east_west(lon_raw)

            data = permute_west_east(data_mask, longrid)#idxwest, idxeast)
            msk = ClimateTools.permute_west_east(msk, longrid)
        else
            data = data_mask
        end

    end

  elseif isempty(poly) # no polygon clipping

      msk = Array{Float64}(ones((size(data_pointer, 1), size(data_pointer, 2))))
      data_ext = ClimateTools.extractdata(data_pointer, msk, idxtimebeg, idxtimeend)
      # replace_missing!(data_ext)
      data_ext = nomissing(data_ext, NaN)
      # data_ext = convert(data_ext, Float32)

    if rotatedgrid

        # Flip data "west-east"
        data = ClimateTools.permute_west_east(data_ext, longrid)

    else
        data = data_ext
    end

  end

  if rotatedgrid
      longrid .= longrid_flip
      lon_raw .= lon_raw_flip
  end

  # Convert units of optional argument data_units is provided
  if data_units == "Celsius" && (vari == "tas" || vari == "tasmax" || vari == "tasmin") && dataunits == "K"
    data .-= 273.15
    dataunits = "°C"
    varattrib["units"] = "Celsius"
    # @warn "Using Celsius can be problematic for arithmetic operations. Best practice is to keep Kelvin and only convert to Celsius at the end with the overloaded ClimateTools.uconvert function."
  end

  if data_units == "mm" && vari == "pr" && (dataunits == "kg m-2 s-1" || dataunits == "mm s-1")

    factor = timeresolution(ds[timename])
    # factor = pr_timefactor(rez)
    data .*= factor.value
    dataunits = "mm"
    varattrib["standard_name"] = "precipitation"
    varattrib["units"] = "mm"
  end

  # Create AxisArray from variable "data"
  if ndims(data) == 3
    # Convert data to AxisArray
    dataOut = AxisArray(data, Axis{Symbol(lonname)}(lon_raw), Axis{Symbol(latname)}(lat_raw), Axis{:time}(timeV))
  elseif ndims(data) == 4 # this imply a 3D field (height component)
    # Get level vector
    dlev = ds[levname][:]
    # Convert data to AxisArray
    dataOut = AxisArray(data, Axis{Symbol(lonname)}(lon_raw), Axis{Symbol(latname)}(lat_raw), Axis{Symbol(levname)}(dlev), Axis{:time}(timeV))
  else
    throw(error("load takes only 3D and 4D variables for the moment"))
  end

  C = ClimGrid(dataOut, longrid=longrid, latgrid=latgrid, msk=msk, grid_mapping=map_attrib, dimension_dict=dimension_dict, timeattrib=timeattrib, model=model, frequency=frequency, experiment=experiment, run=runsim, project=project, institute=institute, filename=file, dataunits=dataunits, latunits=latunits, lonunits=lonunits, variable=vari, typeofvar=vari, typeofcal=caltype, varattribs=varattrib, globalattribs=attribs)

  close(ds)
  # NetCDF.close(ncfile)

  return C


end

"""
    load(files::Array{String,1}, vari::String; poly = ([]), start_date::Date = Date(-4000), end_date::Date = Date(-4000), data_units::String = "")

Loads and merge the files contained in the arrar files.
"""
function load(files::Array{String,1}, vari::String; poly=([]), start_date::Tuple=(Inf,), end_date::Tuple=(Inf,), data_units::String="")

    nfiles = length(files)
    C = Array{ClimGrid}(undef, nfiles) # initialize # TODO better initialization
    datesort = Array{Any}(undef, nfiles)
    Cout = []

    p = Progress(nfiles*2, 3, "Loading files: ")

    for ifile = 1:nfiles

        C[ifile] = load(files[ifile], vari, poly = poly, start_date=start_date, end_date=end_date, data_units=data_units)
        datesort[ifile] = get_timevec(C[ifile])[1]

        next!(p)
    end

    # Sort files based on timevector to ensure that merge results in amonotone increase in time
    idx = sortperm(datesort)
    C = C[idx]

    for imod = 1:nfiles
        if imod == 1
            Cout = C[imod]
        else
            Cout = merge(Cout, C[imod])
        end
        next!(p)
    end

    return Cout
end


"""
    load2D(file::String, vari::String; poly=[], data_units::String="")

Returns a 2D array. Should be used for *fixed* data, such as orography.
"""
function load2D(file::String, vari::String; poly=[], data_units::String="")
    # Get attributes

    ds = NCDatasets.Dataset(file)
    attribs_dataset = ds.attrib
    attribs = Dict(attribs_dataset)

    # Data pointer
    data_pointer = ds[vari]

    project = ClimateTools.project_id(attribs_dataset)
    institute = ClimateTools.institute_id(attribs_dataset)
    model = ClimateTools.model_id(attribs_dataset)
    experiment = ClimateTools.experiment_id(attribs_dataset)
    frequency = ClimateTools.frequency_var(attribs_dataset)
    runsim = ClimateTools.runsim_id(attribs_dataset)
    grid_mapping = ClimateTools.get_mapping(keys(ds))

    if grid_mapping == "Regular_longitude_latitude"
        latstatus = false
    else
        latstatus = true
    end

    # Get dimensions names
    lonname = get_dimname(ds, "X")
    latname = get_dimname(ds, "Y")

    dataunits = ds[vari].attrib["units"]
    latunits = ds[latname].attrib["units"]
    lonunits = ds[lonname].attrib["units"]

    # Create dict with latname and lonname
    dimension_dict = Dict(["lon" => lonname, "lat" => latname])

    # lat_raw = NetCDF.ncread(file, latname)
    lat_raw = nomissing(ds[latname][:], NaN)
    # lon_raw = NetCDF.ncread(file, lonname)
    lon_raw = nomissing(ds[lonname][:], NaN)

    # Get variable attributes
    varattrib = Dict(ds[vari].attrib)

    if latstatus # means we don't have a "regular" grid
        # Get names of grid
        latgrid_name = latgridname(ds)
        longrid_name = longridname(ds)

        # latgrid = NetCDF.ncread(file, latgrid_name)
        latgrid = nomissing(ds[latgrid_name][:], NaN)
        # longrid = NetCDF.ncread(file, longrid_name)
        longrid = nomissing(ds[longrid_name][:], NaN)

        map_attrib = build_grid_mapping(ds, grid_mapping)
        varattrib["grid_mapping"] = grid_mapping

    else # if no grid provided, create one
        longrid, latgrid = ndgrid(lon_raw, lat_raw)
        map_attrib = Dict(["grid_mapping" => grid_mapping])
        varattrib["grid_mapping"] = grid_mapping
    end

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
      data_ext = ClimateTools.extractdata2D(data_pointer, msk)
      # replace_missing!(data_ext)
      data_ext = nomissing(data_ext, NaN)
      # data_ext = convert(data_ext, Float32)

      begin
        I = Base.findall(!isnan, msk)
        idlon, idlat = (getindex.(I, 1), getindex.(I, 2))
      end

      minXgrid = minimum(idlon)
      maxXgrid = maximum(idlon)
      minYgrid = minimum(idlat)
      maxYgrid = maximum(idlat)
      msk = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]
      data_mask = applymask(data_ext, msk) # needed when polygon is not rectangular

      if rotatedgrid
          # Regrid msk to shifted grid if the grid have been rotated to get final mask
          #shiftarray_west_east(msk, longrid)
      end

      # Get lon_raw and lat_raw for such region
      lon_raw = lon_raw[minXgrid:maxXgrid]
      lat_raw = lat_raw[minYgrid:maxYgrid]

      if map_attrib["grid_mapping"] == "Regular_longitude_latitude"

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

              data = permute_west_east(data_mask, longrid)#idxwest, idxeast)
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

          else
            data = data_mask
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

              data = permute_west_east(data_mask, longrid)#idxwest, idxeast)
              msk = ClimateTools.permute_west_east(msk, longrid)

          else
            data = data_mask
          end
      end

    elseif isempty(poly) # no polygon clipping
        msk = Array{Float64}(ones((size(data_pointer, 1), size(data_pointer, 2))))
        data_ext = extractdata2D(data_pointer, msk)
        # replace_missing!(data_ext)
        data_ext = nomissing(data_ext, NaN)
        # data_ext = convert(data_ext, Float32)

      if rotatedgrid
          # Flip data "west-east"
          data = permute_west_east(data_ext, longrid)
      else
        data = data_ext
      end

    end

    if rotatedgrid
        longrid = longrid_flip
        lon_raw = lon_raw_flip
    end

    # Convert data to AxisArray
    dataOut = AxisArray(data, Axis{Symbol(lonname)}(lon_raw), Axis{Symbol(latname)}(lat_raw))



    C = ClimGrid(dataOut, longrid=longrid, latgrid=latgrid, msk=msk, grid_mapping=map_attrib, dimension_dict=dimension_dict, model=model, frequency=frequency, experiment=experiment, run=runsim, project=project, institute=institute, filename=file, dataunits=dataunits, latunits=latunits, lonunits=lonunits, variable=vari, typeofvar=vari, typeofcal="fixed", varattribs=varattrib, globalattribs=attribs)

    close(ds)

    return C

end

model_id(attrib::NCDatasets.Attributes) = get(attrib,"model_id", get(attrib, "parent_source_id", get(attrib,"model","N/A")))

experiment_id(attrib::NCDatasets.Attributes) = get(attrib,"experiment_id",get(attrib,"experiment","N/A"))

project_id(attrib::NCDatasets.Attributes) = get(attrib,"project_id", get(attrib, "mip_era", get(attrib,"project","N/A")))

institute_id(attrib::NCDatasets.Attributes) = get(attrib,"institute_id",get(attrib, "institution_id", get(attrib,"institute","N/A")))

frequency_var(attrib::NCDatasets.Attributes) = get(attrib,"frequency","N/A")

runsim_id(attrib::NCDatasets.Attributes) = get(attrib, "parent_experiment_rip", get(attrib,"driving_model_ensemble_member","N/A"))


"""
    getdim_lat(ds::NCDatasets.Dataset)

Returns the name of the "latitude" dimension and the status related to a regular grid. The latitude dimension is usually "latitude", "lat", "y", "yc", "rlat".
"""
function getdim_lat(ds::NCDatasets.Dataset)

    if sum(keys(ds.dim) .== "rlat") == 1
        return "rlat", true
    elseif sum(keys(ds.dim) .== "lat") == 1
        return "lat", false
    elseif sum(keys(ds.dim) .== "latitude") == 1
        return "latitude", false
    elseif sum(keys(ds.dim) .== "y") == 1
        return "y", true
    elseif sum(keys(ds.dim) .== "yc") == 1
        return "yc", true
    elseif sum(keys(ds.dim) .== "lat_c") == 1
        return "lat_c", false
    else
        error("Manually verify x/lat dimension name")
    end

end

"""
    getdim_lon(ds::NCDatasets.Dataset)

Returns the name of the "longitude" dimension and the status related to a regular grid. The longitude dimension is usually "longitue", "lon", "x", "xc", "rlon".
"""
function getdim_lon(ds::NCDatasets.Dataset)

    if sum(keys(ds.dim) .== "rlon") == 1
        return "rlon", true
    elseif sum(keys(ds.dim) .== "lon") == 1
        return "lon", false
    elseif sum(keys(ds.dim) .== "longitude") == 1
        return "longitude", false
    elseif sum(keys(ds.dim) .== "x") == 1
        return "x", false
    elseif sum(keys(ds.dim) .== "xc") == 1
        return "xc", false
    elseif sum(keys(ds.dim) .== "lon_c") == 1
        return "lon_c", false
    else
        error("Manually verify x/lat dimension name")
    end

end

"""
    latgridname(ds::NCDatasets.Dataset)

Returns the name of the latitude grid when datasets is not on a rectangular grid.
"""
function latgridname(ds::NCDatasets.Dataset)

    # if in("lat", keys(ds))
    #     return "lat"
    # elseif in("latitude", keys(ds))
    #     return "latitude"
    # elseif in("lat_c", keys(ds))
    #     return "lat_c"
    # else
    #     error("Variable name is not supported. File an issue on https://github.com/Balinus/ClimateTools.jl/issues")
    # end

    names = ["degree_north",
             "degrees_north",
             "degreeN",
             "degreesN",
             "degree_N",
             "degrees_N"]

    varnames = keys(ds)

    found_var = "NA"

    for ivar in varnames

        if isdefined(ds[ivar], :attrib)
            if haskey(ds[ivar].attrib, "units")
                if in(ds[ivar].attrib["units"], names)
                    found_var = ivar
                end
            end
        end

    end

    return found_var

end

"""
    longridname(ds::NCDatasets.Dataset)

Returns the name of the longitude grid when datasets is not on a rectangular grid.
"""
function longridname(ds::NCDatasets.Dataset)

    # if in("lon", keys(ds))
    #     return "lon"
    # elseif in("longitude", keys(ds))
    #     return "longitude"
    # elseif in("lon_c", keys(ds))
    #     return "lon_c"
    # else
    #     error("Variable name is not supported. File an issue on https://github.com/Balinus/ClimateTools.jl/issues")
    # end

    names = ["degree_east",
             "degrees_east",
             "degreeE",
             "degreesE",
             "degree_E",
             "degrees_E"]

    varnames = keys(ds)

    found_var = "NA"

    for ivar in varnames

        if isdefined(ds[ivar], :attrib)
            if haskey(ds[ivar].attrib, "units")
                if in(ds[ivar].attrib["units"], names)
                    found_var = ivar
                end
            end
        end

    end

    return found_var



end

"""
    getdim_tim(ds::NCDatasets.Dataset)

Returns the name of the "time" dimension.
"""
function getdim_tim(ds::NCDatasets.Dataset)

    if sum(keys(ds.dim) .== "time") == 1
        return "time", false
    elseif sum(keys(ds.dim) .== "tim") == 1
        return "tim", false
    else
        error("Manually verify time dimension name")
    end

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
        dataout = data[minXgrid:maxXgrid, minYgrid:maxYgrid, idxtimebeg:idxtimeend]
        # Permute dims
        # data = permutedims(data, [3, 1, 2])
    elseif ndims(data) == 4
        dataout = data[minXgrid:maxXgrid, minYgrid:maxYgrid, :, idxtimebeg:idxtimeend]
        # Permute dims
        # data = permutedims(data, [4, 1, 2, 3])
    end

    return dataout
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


"""
    get_mapping(ds::Array{String,1})

Returns the grid_mapping of Dataset *ds*
"""
function get_mapping(K::Array{String,1})

    if in("rotated_pole", K)
        return "rotated_pole"

    elseif in("lambert_conformal_conic", K)
        return "lambert_conformal_conic"

    elseif in("rotated_latitude_longitude", K)
        return "rotated_latitude_longitude"

    elseif in("rotated_mercator", K)
        return "rotated_mercator"

    elseif in("crs", K)
        return "crs"

    elseif in("polar_stereographic", K)
        return "polar_stereographic"

    else
        return "Regular_longitude_latitude"
    end

end

"""
    get_mapping(ds::Array{String,1})

Returns the grid_mapping of Dataset *ds*
"""
function get_mapping(ds::NCDatasets.Dataset, vari)

    grid = ""

    try
        grid = ds[vari].attrib["grid_mapping"]
    catch
        grid = "Regular_longitude_latitude"
    end

    return grid


end

function build_grid_mapping(ds::NCDatasets.Dataset, grid_mapping::String)

    if ClimateTools.@isdefined grid_mapping
        # map_dim = varattrib[grid_mapping]
        map_attrib = Dict(ds[grid_mapping].attrib)
        map_attrib["grid_mapping"] = grid_mapping
    else
        error("File an issue on https://github.com/Balinus/ClimateTools.jl/issues to get the grid supported")
    end

    return map_attrib

end
