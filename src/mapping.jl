"""
    mapclimgrid(C::ClimGrid; region::String="auto", poly, level, mask, caxis, start_date::Tuple, end_date::Tuple, titlestr::String, surface::Symbol, cm::String="", ncolors::Int, center_cs::Bool, filename::String, cs_label::String)

Maps the time-mean average of ClimGrid C. If a filename is provided, the figure is saved in a png format.

Optional keyworkd includes precribed regions (keyword *region*, see list below), spatial clipping by polygon (keyword *poly*) or mask (keyword *mask*, an array of NaNs and 1.0 of the same dimension as the data in ClimGrid C), start_date and end_date. For 4D data, keyword *level* is used to map a given level (defaults to 1). *caxis* is used to limit the colorscale. *cm* is used to manually set the colorscale (see Python documentation for native colorscale keyword), *ncolors* is used to set the number of color classes (defaults to 12). Set *center_cs* to true to center the colorscale (useful for divergent results, such as anomalies, positive/negative temprature). *cs_label* is used for custom colorscale label.

## Arguments for keyword *region* (and shortcuts)
- Europe ("EU")
- NorthAmerica ("NA")
- Canada ("CA")
- Quebec, QuebecNSP ("QC", "QCNSP")
- Americas ("Ams")
- World, WorldAz, WorldEck4 ("W", "Waz", "Weck4")
- Greenwich ("Gr")

## Arguments for keyword *surface*
- :contour
- :contourf
- :pcolormesh
"""
function mapclimgrid(C::ClimGrid; region::String="auto", states::Bool=false, poly=[], level=1, mask=[], caxis=[], start_date::Tuple=(Inf,), end_date::Tuple=(Inf,), titlestr::String="", surface::Symbol=:contourf, cm::String="", ncolors::Int=12, center_cs::Bool=false, filename::String="", cs_label::String="")

  # TODO Add options for custom region

  # Some tests
  if !isempty(mask)
      @assert (size(C[1], 1),size(C[1], 2))==(size(mask, 1),size(mask, 2))
  end

  # get boundaries and lat-lon vectors
  llon = minimum(C.longrid)
  rlon = maximum(C.longrid)
  slat = minimum(C.latgrid)
  nlat = maximum(C.latgrid)

  # Time limits
  if ndims(C[1]) > 2
      timeV = get_timevec(C)
      f = typeof(timeV[1])
      timebeg, timeend = ClimateTools.timeindex(timeV, start_date, end_date, C.frequency, f)
  end

  # =============
  # Colorscale
  # TODO replace C[10] comparison with C.varattribs["standard_name"]
  if isempty(cm)
    if C[10] == "pr" || C[10]=="huss"
      # cm = "YlGnBu"
      cm = cmocean[:cm][:deep]
    elseif C[10]=="tasmax" || C[10]=="tasmin" || C[10]=="tas" || C[10]=="tmax" || C[10]=="tmin" || C[10] == "wbgtmean" || C[10] == "wbgtmax"
      cm = "RdYlBu_r"
    elseif C[10]=="psl" || C[10]=="vp" # pressure
      cm = cmocean[:cm][:deep_r]
    elseif C[10]=="ua" # wind
      cm = cmocean[:cm][:balance]
    else
      cm = "viridis"
    end
  end

  # overide colorscale if the user specify to center the colorscale
  if center_cs
      cm = "RdBu_r"
  end

  if surface == :contourf || surface == :contour
      N = 256
  elseif surface == :pcolormesh
      N = ncolors
  end
  cmap = mpl[:cm][:get_cmap](cm)
  colorlist = cmap(range(0, stop=1, length=N))
  # #
  cm = mpl[:colors][:LinearSegmentedColormap][:from_list]("cm_custom", colorlist, N)

  # =================
  # PLOT DATA
  # Time-average
  if ndims(C[1]) > 2
      data2 = ClimateTools.timeavg(C, timebeg, timeend, mask, poly, level)
  elseif ndims(C[1]) == 2
      data2 = C[1].data
  end

  # Get colorscale limits
  vmin, vmax = ClimateTools.getcslimits(caxis, data2, center_cs)

  # Empty-map generator
  status, fig, ax, m = mapclimgrid(region=region, states=states, llon=llon, rlon=rlon, slat=slat, nlat=nlat)

  x, y = m(C.longrid, C.latgrid) # convert longrid and latgrid to projected coordinates
  if surface == :contourf
    # try
        cs = m[surface](x, y, data2, ncolors, cmap = cm, vmin=vmin, vmax=vmax)
    # catch
    #     cs = m[surface](x, y, data2, ncolors, cmap = cm, vmin=vmin, vmax=vmax)
    # end
  else
    cs = m[surface](x, y, data2, cmap = cm, vmin=vmin, vmax=vmax)
  end

  # Colorbar
  if isempty(cs_label)
      cs_label = getunitslabel(C)
  end
  cbar = colorbar(cs, orientation = "vertical", shrink = 0.7, label=cs_label)

  # Title
  if isempty(titlestr)
      titlestr = titledef(C)
  end
  title(titlestr)

  # Save to "filename" if not empty
  if !isempty(filename)
      PyPlot.savefig(filename, dpi=300)
  end

  return true, fig, ax, cbar
end

"""
    mapclimgrid(; region::String="auto", poly, level, mask, caxis, start_date::Date, end_date::Date)

Empty map generator, when called without a ClimGrid as the positional argument.
"""
function mapclimgrid(;region::String="auto", states::Bool=true, llon=[], rlon=[], slat=[], nlat=[])

    fig, ax = subplots()

    if lowercase(region) == "canada" || lowercase(region) == "ca"
        m = basemap[:Basemap](projection = "lcc", resolution = "l", width=6500000,height=5000000, lat_0 = 62, lon_0 = -95, lat_1 = 45., lat_2 = 55, rsphere = (6378137.00, 6356752.3142))

    elseif lowercase(region) == "quebec" || lowercase(region) == "qc"
        m = basemap[:Basemap](projection = "lcc", resolution = "l", llcrnrlon = -80.5, llcrnrlat = 41., urcrnrlon = -50.566, urcrnrlat = 62.352, lon_0 = -70, lat_1 = 50, rsphere = (6378137.00, 6356752.3142))

    elseif lowercase(region) == "south_quebec" || lowercase(region) == "south_qc"
        m = basemap[:Basemap](projection = "lcc", resolution = "l", llcrnrlon = -75.9, llcrnrlat = 42.6, urcrnrlon = -61.5, urcrnrlat = 49.5, lon_0 = -70, lat_1 = 50, rsphere = (6378137.00, 6356752.3142))

    elseif lowercase(region) == "quebecnsp" || lowercase(region) == "qcnsp"
        m = basemap[:Basemap](projection = "nsper", resolution= "l", satellite_height = 2000000, lon_0 = -72.5, lat_0 = 55)

    elseif lowercase(region) == "americas" || lowercase(region) == "ams"
        m = basemap[:Basemap](projection = "omerc", resolution = "c", width=14000000, height=17000000, lon_0 = -100, lat_0 =    15, lon_1 = -45, lon_2 = -120, lat_1 = -55, lat_2 = 70)

    elseif lowercase(region) == "arctic" || lowercase(region) == "aps"
        m = basemap[:Basemap](projection = "npstere", resolution = "l", boundinglat = 47, lon_0 = 255)

    elseif lowercase(region) == "antarctic" || lowercase(region) == "aaps"
        m = basemap[:Basemap](projection = "spstere", resolution = "l", boundinglat = -60, lon_0 = 210)

    elseif lowercase(region) == "greenwich" || lowercase(region) == "gr"
        m = basemap[:Basemap](projection = "omerc", resolution = "c", width=9000000, height=15000000, lon_0 = 10, lat_0 = 25, lon_1 = -10, lon_2 = 20, lat_1 = -75, lat_2 = 30)

    elseif lowercase(region) == "europe" || lowercase(region) == "eu"
        m = basemap[:Basemap](projection = "lcc", resolution = "l", width = 6800000, height = 4700000, lat_0 = 53, lon_0 = 20, lat_1 = 33, lat_2 = 50, rsphere = (6378137.00, 6356752.3142))

    elseif lowercase(region) == "northamerica" || lowercase(region) == "na"
        m = basemap[:Basemap](projection = "lcc", resolution = "l", llcrnrlon = -135.5, llcrnrlat = 1., urcrnrlon = -10.566, urcrnrlat = 46.352, lon_0 = -95, lat_1 = 50, rsphere = (6378137.00, 6356752.3142))

    elseif lowercase(region) == "southamerica" || lowercase(region) == "sa"
        m = basemap[:Basemap](projection = "lcc", resolution = "l", llcrnrlon = -110., llcrnrlat = -55., urcrnrlon = -30., urcrnrlat = 17., lon_0 = -60, lat_1 = -20, rsphere = (6378137.00, 6356752.3142))

    elseif lowercase(region) == "world" || lowercase(region) == "w"
        m = basemap[:Basemap](projection = "cyl", resolution = "c", llcrnrlat = -90, urcrnrlat = 90, llcrnrlon = -180, urcrnrlon = 180)

    elseif lowercase(region) == "worldaz" || lowercase(region) == "waz"
        m = basemap[:Basemap](projection = "aeqd", resolution = "c", width = 28000000, height = 28000000, lon_0 = -75, lat_0 = 45)

    elseif lowercase(region) == "worldeck4" || lowercase(region) == "weck4"
        m = basemap[:Basemap](projection="eck4", resolution = "c", lon_0 = 0)

    elseif lowercase(region) == "outaouais"
        m = basemap[:Basemap](projection="lcc", resolution="h", llcrnrlon=-78.5, llcrnrlat=45., urcrnrlon=-73.866, urcrnrlat=48.0, lon_0=-75, lat_1=44, rsphere=(6378137.00, 6356752.3142))

    elseif region == "auto"
        m = basemap[:Basemap](projection="cyl", resolution = "c", llcrnrlat = slat, urcrnrlat = nlat, llcrnrlon = llon, urcrnrlon = rlon)
    end

    # Draw the map properties
    m[:drawcoastlines](linewidth = 0.6)
    m[:drawcountries](linewidth = 0.6)

    if states
        m[:drawstates](linewidth = 0.2)
    end

    if lowercase(region) != "quebecnsp" || lowercase(region) != "qcnsp"
      if lowercase(region) == "antarctic" || lowercase(region) == "aaps"
        # Is there a Julia eqvlnt to Py's 'if object not in list'?
          m[:drawmeridians](0:30:360.0, labels = [1,1,1,1], fontsize = 8, linewidth = 0)
      else
          m[:drawparallels](-90:10.0:90, labels = [1,0,0,0], fontsize = 8, linewidth = 0.6)
          m[:drawmeridians](0:30:360.0, labels = [0,0,0,1], fontsize = 8, linewidth = 0.5)
      end
    end

    return true, fig, ax, m

end

"""
    plot(C::ClimGrid, titlefig::String, gridfig::Bool, label::String, color, lw, linestyle)

Plots the spatial average timeserie of ClimGrid `C`.
"""
function PyPlot.plot(C::ClimGrid; poly=[], start_date::Tuple=(Inf,), end_date::Tuple=(Inf,), titlestr::String="", gridfig::Bool=true, label::String="", lw=1.5, linestyle="-", xlimit=[], ylimit=[])

    if !isempty(poly)
        C = spatialsubset(C, poly)
    end
    if !isinf(start_date[1]) || !isinf(end_date[1])
        C = temporalsubset(C, start_date, end_date)
    end

    data = C[1].data
    timevec = get_timevec(C)

    if typeof(timevec[1]) != Date
        if typeof(timevec[1]) == Dates.Year
            timevec = Date.(timevec)
        end
    end

    average = fill(NaN, length(timevec))

    # Spatial mean for each timestep
    for t in 1:length(timevec)
        datatmp = data[:, :, t]
        average[t] = Statistics.mean(datatmp[.!isnan.(datatmp)])
    end

    # figh, ax = subplots()

    if isempty(label)
        label = C.model
    end

    # PLOTTING
    figh = plot(timevec, average, lw=lw, label=label, linestyle=linestyle)
    if !isempty(xlimit)
        xlim(xlimit[1], xlimit[2])
    end
    if !isempty(ylimit)
        ylim(ylimit[1], ylimit[2])
    end
    xlabel("Time")
    ylabel(C.dataunits)
    legend()
    # xticks(timevec[1:20:end])
    if isempty(titlestr)
        titlestr = C.variable
    end
    title(titlestr)
    if gridfig
        grid("on")
    end

    return true, figh

end

"""
    getcslimits(caxis, data, C)

Returns minimum and maximum values of the colorscale axis. Used internally by [`mapclimgrid`](@ref).
"""
function getcslimits(caxis, data, center_cs)

    if !isempty(caxis)
        vmin = caxis[1]
        vmax = caxis[2]
    else
        vmin=minimum(data[.!isnan.(data)])
        vmax=maximum(data[.!isnan.(data)])
    end

    if center_cs # Hack, we want to center the divergent colorscale
        vmin = -vmax
    end

    return vmin, vmax


end

"""
    timeavg(C, timebeg, timeend, mask, poly, level)

Returns an array for mapping purpose. Used internally by [`mapclimgrid`](@ref).
"""
function timeavg(C, timebeg, timeend, mask, poly, level)
    data2 = Array{Float64}(undef, size(C[1], 2), size(C[1], 3))
    if ndims(C[1]) == 3
      data2 = Array(dropdims(Statistics.mean(C[1][:, :, timebeg:timeend], dims=3), dims=3)) #time mean

      # TODO throw error/warning if no grid point inside polygon or mask

      if !isempty(mask) # if mask is already provided
          data2 = data2 .* mask
      elseif !isempty(poly) # if mask needs to be calculated from polygon
          msk = inpolygrid(C.longrid, C.latgrid, poly)
          data2 = data2 .* msk
      end

    # 4D fields
  elseif ndims(C[1]) == 4 # 4D field
      data2 = Array(dropdims(Statistics.mean(C[1][:, :, level, timebeg:timeend], dims=4), dims=3)) # time mean over "level"

      if !isempty(poly)
          msk = inpolygrid(C.longrid, C.latgrid, poly)
          data2 = data2 .* msk
      elseif !isempty(mask)
          data2 = data2 .* mask
      end

    end
    return data2
end

"""
    titledef(C::ClimGrid)

Returns the title. Used internally by [`mapclimgrid`](@ref).
"""
function titledef(C::ClimGrid)
    if ndims(C[1]) > 2

        if typeof((C[1][Axis{:time}][1])) <: AbstractCFDateTime
            begYear = string(Dates.year(C[1][Axis{:time}][1]))
            endYear = string(Dates.year(C[1][Axis{:time}][end]))
        elseif typeof((C[1][Axis{:time}][1])) == Dates.Year
            begYear = string(C[1][Axis{:time}][1])[1:4]
            endYear = string(C[1][Axis{:time}][end])[1:4]
        elseif typeof((C[1][Axis{:time}][1])) == Int
            begYear = string(C[1][Axis{:time}][1])[1:4]
            endYear = string(C[1][Axis{:time}][end])[1:4]
        end
        titlestr = string(C[3], " - ", C[5], " - ", C.variable, " - ", begYear, " - ", endYear)

    else
        titlestr = string(C[3], " - ", C[5], " - ", C.variable)
    end

    return titlestr

end

"""
    getunitslabel(C::ClimGrid)

Return verbose label for colorbar. Used internally by [`mapclimgrid`](@ref).
"""
function getunitslabel(C::ClimGrid)

    try
        standardname = C.varattribs["standard_name"]
    catch
        standardname = C.varattribs["long_name"]
    end

    units = Dict(["air_temperature" => "Air temperature",
    "specific_humidity" => "Specific humidity (%)",
    "precipitation" => "Precipitation",
    "precipitation_flux" => "Precipitation flux",
    "air_pressure_at_sea_level" => "Air pressure at sea level",
    "eastward_wind" => "Eastward wind",
    "daily maximum temperature" => "Air temperature"
    ])

    label = ""

    try
        label = string(units[standardname], " (", C[2], ")")
    catch
        label = C[2]
    end

    return label

end
