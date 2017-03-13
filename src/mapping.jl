function mapclimgrid(C::ClimGrid; region::String = "auto")

  # get boundaries and lat-lon vectors

  llon = minimum(C.data[Axis{:lon}][:])
  rlon = maximum(C.data[Axis{:lon}][:])
  slat = minimum(C.data[Axis{:lat}][:])
  nlat = maximum(C.data[Axis{:lat}][:])
  lat = C[1][Axis{:lat}][:]
  lon = C[1][Axis{:lon}][:]


  if (sum( (C.data[Axis{:lon}][:] .> 355) & (C.data[Axis{:lon}][:] .< 5)) > 0)
    if (sum( (C.data[Axis{:lon}][:] .< 185) & (C.data[Axis{:lon}][:] .> 175) ) == 0)
      (llon, rlon) = (rlon, llon)
    end
  end
  figure(figsize=(12,7), dpi = 120)
  if region == "Canada"
    m = basemap[:Basemap](width=6500000,height=5000000, rsphere = (6378137.00, 6356752.3142),      resolution = "l", area_thresh = 800., projection = "lcc", lat_1 = 45., lat_2 = 55, lat_0 = 62, lon_0 = -95.)
  elseif region == "World"
    m = basemap[:Basemap](projection="cyl", llcrnrlat = -90, urcrnrlat = 90, llcrnrlon = 0, urcrnrlon = 360, resolution = "c")
  elseif region == "auto"
    m = basemap[:Basemap](projection="cyl", llcrnrlat = slat, urcrnrlat = nlat, llcrnrlon = llon, urcrnrlon = rlon, resolution = "c")
  end

  # Plot the value of the ClimGrid data
  # m[:drawlsmask](land_color = "#00441b", ocean_color = "#8be5e5", lakes = true)
  m[:drawcoastlines]()
  m[:drawstates]()
  m[:drawcountries]()
  m[:drawmeridians](0:30:360.0, labels = [0,0,0,1], fontsize = 10);
  m[:drawparallels](-90:10.0:90, labels = [1,0,0,0], fontsize = 10);
  # parallels = np.arange(llon, rlon, 20.)
  # m[:drawparallels](parallels, labels = [1,0,0,0], fontsize = 8)
  # meridians = np.arange(slat, nlat, 20.)
  # m[:drawmeridians](meridians, labels = [0,0,0,1], fontsize = 8)

  lon2, lat2 = np[:meshgrid](lon, lat)
  x, y = m(lon2, lat2)

  if C[10] == "pr"
    cm = "YlGnBu"
  elseif C[10] == "tasmax" || C[10] == "tasmin" || C[10] == "tas"
    cm = "YlOrBr"
  else
    cm = "viridis"
  end

  if length(size(C[1])) > 2

    cs = m[:contourf](x, y, squeeze(mean(convert(Array, C[1]),1),1)', cmap=get_cmap(cm))
  else
    cs = m[:contourf](x, y, convert(Array,C[1][:,:])', cmap=get_cmap(cm))
  end

  cbar = colorbar(cs, orientation = "vertical", shrink = 0.5, label = C[2])
  # cbar[:set_label] = C.dataunits
  title(string(C[3], "-", C[4], "-", C[5], " - ", C.var))

  return 1
end

function drawmap(m, C::ClimGrid; NoLands::Bool = false)
  if (NoLands)
    m[:drawcoastlines](color="#555555")
  else
    m[:fillcontinents](color="#555555")
  end
  m[:drawmeridians](0:30:360.0, labels = [0,0,0,1], fontsize = 10);
  m[:drawparallels](-90:10.0:90, labels = [1,0,0,0], fontsize = 10);
end

function lonrotate(plon :: Array{Float64,2}, rr)
  if ((plon[1,1] > 180) & (plon[end,1] < 180))
    lon = mean(plon ,2)
    return [rr[lon .<= 180,:] ; rr[lon .> 180,:]]
  else
    return rr
  end
end
