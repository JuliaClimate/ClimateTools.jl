function mapclimgrid(C::ClimGrid; region::String = "auto", poly = [], level = 1, mask = [])

  # TODO Add options for custom time period as input and custom region

  # Some tests
  if !isempty(mask)
      @assert (size(C[1], 2), size(C[1], 3)) == (size(mask, 1), size(mask, 2))
      mask = mask'
  end

  # get boundaries and lat-lon vectors
  llon = minimum(C.data[Axis{:lon}][:])
  rlon = maximum(C.data[Axis{:lon}][:])
  slat = minimum(C.data[Axis{:lat}][:])
  nlat = maximum(C.data[Axis{:lat}][:])
  lat = C[1][Axis{:lat}][:]
  lon = C[1][Axis{:lon}][:]

  # if the data is from a GCM, we sometimes needs to extend the lat-lon and data to avoid "white space"
  if rlon > 355 && llon < 5 && rlon < 359.99
    rlon = 360; llon = 0;
    push!(lon, 359.99)
  end

  if (sum( (C.data[Axis{:lon}][:] .> 355) .& (C.data[Axis{:lon}][:] .< 5)) > 0)
    if (sum( (C.data[Axis{:lon}][:] .< 185) .& (C.data[Axis{:lon}][:] .> 175) ) == 0)
      (llon, rlon) = (rlon, llon)
    end
  end
  # figh = figure(figsize=(12,7), dpi = 120)
  figh, ax = subplots(figsize=(6, 4))
  if region == "Canada"
      m = basemap[:Basemap](width=6500000,height=5000000, rsphere = (6378137.00, 6356752.3142),      resolution = "c", projection = "lcc", lat_1 = 45., lat_2 = 55, lat_0 = 62, lon_0 = -95.)

  elseif region == "Quebec"
      m = basemap[:Basemap](llcrnrlon = -80.5, llcrnrlat = 41., urcrnrlon = -50.566, urcrnrlat = 62.352, rsphere = (6378137.00, 6356752.3142), resolution = "c", projection = "lcc",  lat_1 = 50., lon_0 = -70.)

  elseif region == "World"
      m = basemap[:Basemap](projection="cyl", llcrnrlat = -90, urcrnrlat = 90, llcrnrlon = 0, urcrnrlon = 360, resolution = "c")

  elseif region == "Europe"
      m = basemap[:Basemap](width=6800000, height = 4500000, rsphere = (6378137.00, 6356752.3142), resolution = "c", projection = "lcc", lat_1 = 30., lat_2 = 45, lat_0 = 52, lon_0 = 10.)

  elseif region == "NorthAmerica"
      m = basemap[:Basemap](llcrnrlon = -135.5, llcrnrlat = 1., urcrnrlon = -10.566, urcrnrlat = 46.352, rsphere = (6378137.00, 6356752.3142), resolution = "c", projection = "lcc",  lat_1 = 50., lon_0 = -107.)

  elseif region == "auto"
      m = basemap[:Basemap](projection="cyl", llcrnrlat = slat, urcrnrlat = nlat, llcrnrlon = llon, urcrnrlon = rlon, resolution = "c")
  end

  # Draw the map properties
  m[:drawcoastlines]()
  m[:drawstates]()
  m[:drawcountries]()
  m[:drawmeridians](0:30:360.0, labels = [0,0,0,1], fontsize = 8);
  m[:drawparallels](-90:10.0:90, labels = [1,0,0,0], fontsize = 8);


  lon2, lat2 = np[:meshgrid](lon, lat)
  x, y = m(lon2, lat2)

  if C[10] == "pr"
    # cm = "YlGnBu"
    cm = cmocean[:cm][:deep]
  elseif C[10] == "tasmax" || C[10] == "tasmin" || C[10] == "tas"
    cm = "YlOrBr"
    # cm = cmocean[:cm][:solar]
  else
    cm = "viridis"
  end

  # -----------------------
  # Plot the data
  # 3D fields
  if ndims(C[1]) == 3
    data = squeeze(mean(convert(Array, C[1]),1),1)' #time mean

    if rlon > 355 && llon < 5 # to avoid white space along the longitude 0 (this is purely cosmetic for maps and should not be used for calculation!)
      data2 = Array{Float64}(size(data, 1), size(data, 2) + 1)
      data2[:, 1:end-1] = data
      data2[:, end] = data[:, 1]

      if !isempty(mask)
          mask2 = Array{Float64}(size(mask, 1), size(mask, 2) + 1)
          mask2[:, 1:end-1] = mask
          mask2[:, end] = mask[:, 1]
          mask = mask2
      end

    else
      data2 = data
    end
    # plot data

    if !isempty(mask) # if mask is already provided
        cs = m[:contourf](x .* mask, y .* mask, data2 .* mask, cmap = get_cmap(cm))
    elseif !isempty(poly) # if mask needs to be calculated from polygon
        msk = inpolyvec(lon, lat, poly)'
        cs = m[:contourf](x .* msk, y .* msk, data2 .* msk, cmap = get_cmap(cm))
    else
        cs = m[:contourf](x, y, data2, cmap = get_cmap(cm))

    end

  # 2D fields
  elseif ndims(C[1]) == 2
      if isempty(poly)
          cs = m[:contourf](x, y, convert(Array,C[1][:,:])', cmap = get_cmap(cm))
      else
          msk = inpolyvec(lon, lat, poly)'
          cs = m[:contourf](x .* msk, y .* msk, convert(Array,C[1][:,:])' .* msk, cmap = get_cmap(cm))
      end

  # 4D fields
  elseif ndims(C[1]) == 4 # 3D field
    data = squeeze(mean(convert(Array, C[1][:, :, :, level]),1),1)' # time mean over "level"
    # data = datatmp[:,:, level]';
    if rlon > 355 && llon < 5
      data2 = Array{Float64}(size(data, 1), size(data, 2) + 1)
      data2[:, 1:end-1] = data
      data2[:, end] = data[:, 1]
      # push!(data, data[:, 1])
      if !isempty(mask)
          mask2 = Array{Float64}(size(mask, 1), size(mask, 2) + 1)
          mask2[:, 1:end-1] = mask
          mask2[:, end] = mask[:, 1]
          mask = mask2
      end
    else
      data2 = squeeze(mean(convert(Array, C[1][:, level, :, :]),1),1)' #time mean over "level"
    end

    if !isempty(poly)
        msk = inpolyvec(lon, lat, poly)'
        cs = m[:contourf](x .* msk, y .* msk, data2 .* msk, cmap=get_cmap(cm))
    elseif !isempty(mask)
        cs = m[:contourf](x .* mask, y .* mask, data2 .* mask, cmap = get_cmap(cm))
    else
        cs = m[:contourf](x, y, data2, cmap = get_cmap(cm))
    end

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

  return true, figh, ax, cbar
end



# function drawmap(m, C::ClimGrid; NoLands::Bool = false)
#   if (NoLands)
#     m[:drawcoastlines](color="#555555")
#   else
#     m[:fillcontinents](color="#555555")
#   end
#   m[:drawmeridians](0:30:360.0, labels = [0,0,0,1], fontsize = 10);
#   m[:drawparallels](-90:10.0:90, labels = [1,0,0,0], fontsize = 10);
# end
#
# function lonrotate(plon :: Array{Float64,1}, rr)
#   if ((plon[1] > 180) & (plon[end] < 180))
#     lon = mean(plon ,2)
#     return [rr[lon .<= 180,:] ; rr[lon .> 180,:]]
#   else
#     return rr
#   end
# end
