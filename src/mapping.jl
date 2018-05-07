function mapclimgrid(C::ClimGrid; region::String = "auto", states::Bool = false, poly = [], level = 1, mask = [], caxis=[], start_date::Date=Date(-4000), end_date::Date=Date(-4000))

  # TODO Add options for custom time period as input and custom region

  # Some tests
  if !isempty(mask)
      @assert (size(C[1], 2),size(C[1], 3))==(size(mask, 1),size(mask, 2))
  end

  # get boundaries and lat-lon vectors
  llon = minimum(C.longrid)
  rlon = maximum(C.longrid)
  slat = minimum(C.latgrid)
  nlat = maximum(C.latgrid)

  # Time limits
  timeV = C[1][Axis{:time}][:]
  timebeg, timeend = timeindex(timeV, start_date, end_date, C.frequency)

  # =============
  # MAP DEFINITION

  # Projection
  fig, ax = subplots(figsize=(8, 6))
  if region == "Canada"
      m = basemap[:Basemap](projection = "lcc", resolution = "l", width=6500000,height=5000000, lat_0 = 62, lon_0 = -95, lat_1 = 45., lat_2 = 55, rsphere = (6378137.00, 6356752.3142))

  elseif region == "Quebec"
      m = basemap[:Basemap](projection = "lcc", resolution = "l", llcrnrlon = -80.5, llcrnrlat = 41., urcrnrlon = -50.566, urcrnrlat = 62.352, lon_0 = -70, lat_1 = 50, rsphere = (6378137.00, 6356752.3142))

  elseif region == "QuebecNSP"
      m = basemap[:Basemap](projection = "nsper", resolution= "l", satellite_height = 2000000, lon_0 = -72.5, lat_0 = 55)

  elseif region == "Americas"
      m = basemap[:Basemap](projection = "omerc", resolution = "c", width=14000000, height=17000000, lon_0 = -100, lat_0 =    15, lon_1 = -45, lon_2 = -120, lat_1 = -55, lat_2 = 70)

  elseif region == "Greenwich"
      m = basemap[:Basemap](projection = "omerc", resolution = "c", width=9000000, height=15000000, lon_0 = 10, lat_0 = 25, lon_1 = -10, lon_2 = 20, lat_1 = -75, lat_2 = 30)

  elseif region == "Europe"
      m = basemap[:Basemap](projection = "lcc", resolution = "l", width = 6800000, height = 4700000, lat_0 = 53, lon_0 = 20, lat_1 = 33, lat_2 = 50, rsphere = (6378137.00, 6356752.3142))

  elseif region == "NorthAmerica"
      m = basemap[:Basemap](projection = "lcc", resolution = "l", llcrnrlon = -135.5, llcrnrlat = 1., urcrnrlon = -10.566, urcrnrlat = 46.352, lon_0 = -107, lat_1 = 50, rsphere = (6378137.00, 6356752.3142))

  elseif region == "World"
      m = basemap[:Basemap](projection = "cyl", resolution = "c", llcrnrlat = -90, urcrnrlat = 90, llcrnrlon = -180, urcrnrlon = 180)

  elseif region == "WorldAz"
      m = basemap[:Basemap](projection = "aeqd", resolution = "c", width = 28000000, height = 28000000, lon_0 = -75, lat_0 = 45)

  elseif region == "WorldEck4"
      m = basemap[:Basemap](projection="eck4", resolution = "c", lon_0 = -210)

  elseif region == "auto"
      m = basemap[:Basemap](projection="cyl", resolution = "c", llcrnrlat = slat, urcrnrlat = nlat, llcrnrlon = llon, urcrnrlon = rlon)
  end

  # Draw the map properties
  m[:drawcoastlines](linewidth = 0.6)
  m[:drawcountries](linewidth = 0.6)

  if states
      m[:drawstates](linewidth = 0.2)
  end

  if region != "QuebecNSP"
      m[:drawparallels](-90:10.0:90, labels = [1,0,0,0], fontsize = 8, linewidth = 0.6)
      m[:drawmeridians](0:30:360.0, labels = [0,0,0,1], fontsize = 8, linewidth = 0.5)
  end


  x, y = m(C.longrid, C.latgrid)

  # Colorscale
  # TODO replace C[10] comparison with C.varattribs["standard_name"]

  if C[10] == "pr" || C[10]=="huss"
      # cm = "YlGnBu"
      cm = cmocean[:cm][:deep]
  elseif C[10]=="tasmax" || C[10]=="tasmin" || C[10]=="tas" || C[10]=="tmax" || C[10]=="tmin"
      cm = "YlOrBr"
      # cm = cmocean[:cm][:solar]
  elseif C[10]=="psl"
      cm = cmocean[:cm][:deep_r]
  else
      cm = "viridis"
  end

  # =================
  # PLOT DATA
  # 3D fields
  if ndims(C[1]) == 3
    data2 = squeeze(mean(C[1][timebeg:timeend, :, :], 1), 1) #time mean

    # TODO throw error/warning if no grid point inside polygon or mask

    if !isempty(mask) # if mask is already provided
        data2 = data2 .* mask
    elseif !isempty(poly) # if mask needs to be calculated from polygon
        msk = inpolygrid(C.longrid, C.latgrid, poly)
        data2 = data2 .* msk
    end

    if !isempty(caxis)
        vmin = caxis[1]
        vmax = caxis[2]
    else
        vmin=minimum(data2[.!isnan.(data2)])
        vmax=maximum(data2[.!isnan.(data2)])
    end

    # Plot data
    cs = m[:contourf](x, y, data2, cmap = get_cmap(cm), vmin=vmin, vmax=vmax)

    # TODO create 2D ClimGrid
  # # 2D fields
  # elseif ndims(C[1]) == 2
  #     if isempty(poly)
  #         cs = m[:contourf](x, y, convert(Array,C[1][:,:])', cmap = get_cmap(cm))
  #     else
  #         msk = inpolyvec(lon, lat, poly)'
  #         cs = m[:contourf](x .* msk, y .* msk, convert(Array,C[1][:,:])' .* msk, cmap = get_cmap(cm))
  #     end

  # 4D fields
  elseif ndims(C[1]) == 4 # 3D field
    data2 = squeeze(mean(C[1][timebeg:timeend, :, :, level], 1), 1) # time mean over "level"

    if !isempty(poly)
        msk = inpolygrid(C.longrid, C.latgrid, poly)
        data2 = data2 .* msk
    elseif !isempty(mask)
        data2 = data2 .* mask
    end
    cs = m[:contourf](x, y, data2, cmap = get_cmap(cm))

  end

  # Colorbar
  cbar = colorbar(cs, orientation = "vertical", shrink = 0.7, label = C[2])

  if ndims(C[1]) > 2

      if typeof((C[1][Axis{:time}][1])) == Date
          begYear = string(Base.Dates.year(C[1][Axis{:time}][1]))
          endYear = string(Base.Dates.year(C[1][Axis{:time}][end]))
      elseif typeof((C[1][Axis{:time}][1])) == Base.Dates.Year
          begYear = string(C[1][Axis{:time}][1])[1:4]
          endYear = string(C[1][Axis{:time}][end])[1:4]
      elseif typeof((C[1][Axis{:time}][1])) == Int
          begYear = string(C[1][Axis{:time}][1])[1:4]
          endYear = string(C[1][Axis{:time}][end])[1:4]
      end
      title(string(C[3], " - ", C[4], " - ", C[5], " - ", C.variable, " - ", begYear, " - ", endYear))
  else
      title(string(C[3], " - ", C[4], " - ", C[5], " - ", C.variable))
  end



  return true, fig, ax, cbar
end


function mapclimgrid(;region::String = "auto", states::Bool = true)

    fig, ax = subplots(figsize=(8, 6))

    if region == "Canada"
        m = basemap[:Basemap](projection = "lcc", resolution = "l", width=6500000,height=5000000, lat_0 = 62, lon_0 = -95, lat_1 = 45., lat_2 = 55, rsphere = (6378137.00, 6356752.3142))

    elseif region == "Quebec"
        m = basemap[:Basemap](projection = "lcc", resolution = "l", llcrnrlon = -80.5, llcrnrlat = 41., urcrnrlon = -50.566, urcrnrlat = 62.352, lon_0 = -70, lat_1 = 50, rsphere = (6378137.00, 6356752.3142))

    elseif region == "QuebecNSP"
        m = basemap[:Basemap](projection = "nsper", resolution= "l", satellite_height = 2000000, lon_0 = -72.5, lat_0 = 55)

    elseif region == "Americas"
        m = basemap[:Basemap](projection = "omerc", resolution = "c", width=14000000, height=17000000, lon_0 = -100, lat_0 =    15, lon_1 = -45, lon_2 = -120, lat_1 = -55, lat_2 = 70)

    elseif region == "Greenwich"
        m = basemap[:Basemap](projection = "omerc", resolution = "c", width=9000000, height=15000000, lon_0 = 10, lat_0 = 25, lon_1 = -10, lon_2 = 20, lat_1 = -75, lat_2 = 30)

    elseif region == "Europe"
        m = basemap[:Basemap](projection = "lcc", resolution = "l", width = 6800000, height = 4700000, lat_0 = 53, lon_0 = 20, lat_1 = 33, lat_2 = 50, rsphere = (6378137.00, 6356752.3142))

    elseif region == "NorthAmerica"
        m = basemap[:Basemap](projection = "lcc", resolution = "l", llcrnrlon = -135.5, llcrnrlat = 1., urcrnrlon = -10.566, urcrnrlat = 46.352, lon_0 = -107, lat_1 = 50, rsphere = (6378137.00, 6356752.3142))

    elseif region == "World"
        m = basemap[:Basemap](projection = "cyl", resolution = "c", llcrnrlat = -90, urcrnrlat = 90, llcrnrlon = -180, urcrnrlon = 180)

    elseif region == "WorldAz"
        m = basemap[:Basemap](projection = "aeqd", resolution = "c", width = 28000000, height = 28000000, lon_0 = -75, lat_0 = 45)

    elseif region == "WorldEck4"
        m = basemap[:Basemap](projection="eck4", resolution = "c", lon_0 = -210)

    elseif region == "auto"
        m = basemap[:Basemap](projection="cyl", resolution = "c", llcrnrlat = slat, urcrnrlat = nlat, llcrnrlon = llon, urcrnrlon = rlon)
    end

    # Draw the map properties
    m[:drawcoastlines](linewidth = 0.6)
    m[:drawcountries](linewidth = 0.6)

    if states
        m[:drawstates](linewidth = 0.2)
    end

    if region != "QuebecNSP"
        m[:drawparallels](-90:10.0:90, labels = [1,0,0,0], fontsize = 8, linewidth = 0.6)
        m[:drawmeridians](0:30:360.0, labels = [0,0,0,1], fontsize = 8, linewidth = 0.5)
    end

    return true, fig, ax

end
