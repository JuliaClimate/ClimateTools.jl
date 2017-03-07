function simplemap(C::ClimGrid)

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
  m = basemap.Basemap(projection="cyl", llcrnrlat = slat, urcrnrlat = nlat, llcrnrlon = llon, urcrnrlon = rlon, resolution = "c")

  # Plot the value of the ClimGrid data
  # m[:drawlsmask](land_color = "#00441b", ocean_color = "#8be5e5", lakes = true)
  m[:drawcoastlines]()
  m[:drawstates]()
  m[:drawcountries]()
  # parallels = np.arange(llon, rlon, 20.)
  # m[:drawparallels](parallels, labels = [1,0,0,0], fontsize = 8)
  # meridians = np.arange(slat, nlat, 20.)
  # m[:drawmeridians](meridians, labels = [0,0,0,1], fontsize = 8)

  lon2, lat2 = np.meshgrid(lon, lat)
  x, y = m(lon2, lat2)
  # test = convert(Array,C[1][3,:,:])'
  if length(size(C[1]))>2
    cs = m[:contourf](x, y, squeeze(mean(convert(Array, C[1]),1),1)')
  else
    cs = m[:contourf](x, y, convert(Array,C[1][:,:])')
  end
  cbar = colorbar(cs, orientation = "vertical", shrink = 0.5)
  cbar[:set_label] = C.dataunits


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
#
# function canvas( pa :: AnArea_Regular,xarr :: Array{Float64,2})
#   return map(x->lonrotate(pa.plon,x),(pa.plon, pa.plat, xarr))
# end
#
# function canvas( pa :: AnArea_Regular,xarr :: Array{Bool,2})
#   (plon,plat,xx) = map( x-> lonrotate(plon,x), (pa.plon, pa.plat, xarr))
#   return map(plon[xx], plat[xx])
# end
#
# function canvas( pa :: AnArea_Regular,xarr :: BitArray{2})
#   (plon,plat,xx) = map( x-> lonrotate(plon,x), (pa.plon, pa.plat, xarr))
#   return map(plon[xx], plat[xx])
# end
#
# function canvas(x :: Array{Float64,1}, y :: Array{Float64,1},  xarr :: Array{Float64,2})
#   (px, py) = begin
#     p1 = x .+ zeros(y')
#     p2 = y .+ zeros(x')
#     (p1, p2')
#   end
#   return (px, py,xarr)
# end
#
# function canvas(x :: Array{Float64,1}, y :: Array{Float64,1},  xarr :: Array{Float64})
#   xarr1 = squeeze(xarr, tuple(findin(size(xarr),1)...))
#   @assert ( ndims(xarr1) == 2)
#   return canvas(x, y, xarr1)
# end

abstract LevelMethods

type CenterLevels <: LevelMethods
  nlevels :: Int64
  rangeCalculator :: Function
end

function CenterLevels(; nlevels :: Int64 = 5, center_by :: Function = x-> mean(x))
  rangeCalculator = function ( x :: Array{Float64})
    center = center_by(x)
    r1 = 3*std(x)
    return (center - r1, center + r1)
  end
  return CenterLevels(nlevels, rangeCalculator)
end
export CenterLevels

type PDFlevels <: LevelMethods
  nlevels :: Int64
  low_cutoff :: Float64
  high_cutoff :: Float64
end

function PDFlevels(; nlevels :: Int64 = 5, low_cutoff :: Float64 = 0.15)
  return PDFlevels( nlevels, low_cutoff, 1-low_cutoff)
end

function levels( xarr :: Array{Float64} , m :: PDFlevels)
  x1 = sort(xarr[:])
  lx = length(xarr)
  low = x1[floor(Int64,m.low_cutoff * lx)]
  high = x1[floor(Int64,m.high_cutoff * lx)]
  rr = collect( low:(high-low)/(m.nlevels):high )
  return ((:levels, rr), (:extend, "both"))
end


function levels( xarr :: Array{Float64} , m :: CenterLevels)
  (low, high) = m.rangeCalculator(xarr)
  rr = collect( low:(high-low)/(m.nlevels):high )
  return ((:levels, rr), (:extend, "both"))
end
