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

    return ClimGrid(dataOut, longrid=longrid, latgrid=latgrid, msk=msk, grid_mapping=C.grid_mapping, timeattrib=C.timeattrib, dimension_dict=C.dimension_dict, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)

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
    extractpoly(file::String; n::Int=1)

Returns the n-th polygon contained in *file*.
"""
function extractpoly(file::String; n::Int=1)

    shp = open(file) do fd
        read(fd, Shapefile.Handle)
    end

    poly = shapefile_coords_poly(shp.shapes[n])
    return poly
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

function meridian_check(poly)

    minpoly = minimum(poly[1, .!isnan.(poly[1, :])])
    maxpoly = maximum(poly[1, .!isnan.(poly[1, :])])

    meridian = false
    if sign(minpoly)*sign(maxpoly) == sign(-1.0) # polygon crosses the meridian
       meridian = true
    end

    return meridian

end

function findmindist(p1::Tuple{Float64, Float64}, p2::Tuple{Array{Float64,1},Array{Float64,1}}; R=6372.8)

    dist = Array{Float64}(undef, length(l2[1]), 1)

    for idx = 1:length(l2[1])
        l3 = (l2[1][idx],l2[2][idx])
        dist[idx] = haversine(l1, l3, 6372.8)
    end

    return findmin(dist)

end
