"""
    write(C::ClimGrid, filename::String)

Write to disk ClimGrid C to netCDF file.
"""
function write(C::ClimGrid, filename::String)

    # Test extension
    if extension(filename) .!= ".nc"
        filename = string(filename, ".nc")
    end

    # This creates a new NetCDF file
    # The mode "c" stands for creating a new file (clobber)
    ds = Dataset(filename, "c")

    latsymbol, lonsymbol = ClimateTools.getsymbols(C)
    x, y, timevec = ClimateTools.getdims(C)
    longrid, latgrid = ClimateTools.getgrids(C)

    # DIMENSIONS
    defDim(ds, string(lonsymbol),length(x))
    defDim(ds, string(latsymbol), length(y))
    defDim(ds, "time", length(timevec))

    nclon = defVar(ds, string(lonsymbol), Float64, (string(lonsymbol),))
    nclon.attrib["units"] = C.lonunits#"degrees_east"
    nclon.attrib["long_name"] = "longitude"
    nclon.attrib["standard_name"] = "longitude"
    nclon.attrib["actual_range"] = [minimum(x), maximum(x)]

    nclat = defVar(ds, string(latsymbol), Float64, (string(latsymbol),))
    nclat.attrib["units"] = C.latunits#"degrees_north"
    nclat.attrib["long_name"] = "latitude"
    nclat.attrib["standard_name"] = "latitude"
    nclat.attrib["actual_range"] = [minimum(y), maximum(y)]

    if latsymbol != :lat # a lon-lat grid is needed
        ncrlat = defVar(ds, "lat", Float64, (string(lonsymbol), string(latsymbol)))
        ncrlat.attrib["long_name"] = "latitude"
        ncrlat.attrib["units"] = "degrees"
        ncrlat.attrib["standard_name"] = "grid_latitude"
        ncrlat.attrib["axis"] = "Y"
        ncrlat.attrib["coordinate_defines"] = "point"
        ncrlat.attrib["actual_range"] = [minimum(latgrid), maximum(latgrid)]
    end

    if lonsymbol != :lon
        ncrlon = defVar(ds,"lon", Float64, (string(lonsymbol), string(latsymbol)))
        ncrlon.attrib["long_name"] = "longitude"
        ncrlon.attrib["units"] = "degrees"
        ncrlon.attrib["standard_name"] = "grid_longitude"
        ncrlon.attrib["axis"] = "X"
        ncrlon.attrib["coordinate_defines"] = "point"
        ncrlon.attrib["actual_range"] = [minimum(longrid), maximum(longrid)]
    end

    if ClimateTools.@isdefined E_futur.grid_mapping["grid_mapping_name"]
        ncmapping = defVar(ds, C.grid_mapping["grid_mapping_name"], Char, ())
    elseif ClimateTools.@isdefined E_futur.grid_mapping["grid_mapping"]
        ncmapping = defVar(ds, C.grid_mapping["grid_mapping"], Char, ())
    end
    for iattr in keys(C.grid_mapping)
        ncmapping.attrib[iattr] = C.grid_mapping[iattr]
    end

    # time
    timeout = Array{Float64, 1}(undef, length(timevec))
    for itime = 1:length(timevec)
        timeout[itime] = NCDatasets.timeencode([timevec[itime]], C.timeattrib["units"], C.timeattrib["calendar"])[1]
    end


    # daysfrom = string("days since ", string(timevec[1] - Dates.Day(1))[1:10], " ", string(timevec[1])[12:end])

    nctime = defVar(ds,"time", Float64, ("time",))
    nctime.attrib["long_name"] = "time"
    nctime.attrib["standard_name"] = "time"
    nctime.attrib["axis"] = "T"
    nctime.attrib["calendar"] = C.typeofcal#"365_day"
    nctime.attrib["units"] = C.timeattrib["units"]#daysfrom
    # nctime.attrib["bounds"] = "time_bnds"
    nctime.attrib["coordinate_defines"] = "point"

    # Define the variables contained in ClimGrid C
    v = defVar(ds, C.variable, Float32, (string(lonsymbol),string(latsymbol), "time"))    

    # write attributes
    for iattr in keys(C.varattribs)
        v.attrib[iattr] = C.varattribs[iattr]
    end

    # Define global attributes
    for iattr in keys(C.globalattribs)
        ds.attrib[iattr] = C.globalattribs[iattr]
    end

    v[:] = C[1].data

    # Time vector
    nctime[:] = timeout

    # Dimensions
    nclon[:] = x
    nclat[:] = y

    # Dimensions
    if latsymbol != :lat
        ncrlon[:] = longrid
    end
    if latsymbol != :lat
        ncrlat[:] = latgrid
    end

    # Close file
    close(ds)

end
