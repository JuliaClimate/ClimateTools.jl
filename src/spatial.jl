"""
    moving_fct(xout, xin, fct)

Used as inner function of moving_fct(cube::YAXArray).

See also [`moving_fct(cube::YAXArray)`](@ref).
"""
function moving_fct(xout, xin; fct::Function=mean)
    xout .= fct(skipmissing(xin))
end

"""
    moving_fct(cube::YAXArray; fct=mean, lat_before=1, lat_after=1, lon_before=1, lon_after=1))

Spatial moving function. lat_before/after and lon_before/after is the number of pixel used for the bounding box and fct is the function applied
"""
function moving_fct(cube::YAXArray; fct::Function=mean, latsouth_pixel=1, latnorth_pixel=1, lonwest_pixel=1, loneast_pixel=1)
    
    indims = InDims(MovingWindow("latitude", latsouth_pixel,latnorth_pixel),
        MovingWindow("longitude", lonwest_pixel,loneast_pixel))
    outdims=OutDims()
    mapCube(moving_fct, cube; fct=fct, indims=indims, outdims=outdims)
end




"""
    spatialsubset(cube::YAXArray, poly)

Return a polygon-masked spatial subset of a YAXArray.

The polygon is expected to be in lon/lat coordinates on a -180 to 180
longitude convention. `poly` can be provided as a 2xN matrix or an Nx2
matrix.

"""
function _normalize_polygon(poly)
    ndims(poly) == 2 || error("Polygon must be a 2D array with lon/lat coordinates.")

    if size(poly, 1) == 2
        return Float64.(poly)
    elseif size(poly, 2) == 2
        return Float64.(poly')
    else
        error("Polygon must be a 2xN or Nx2 array.")
    end
end

function _find_spatial_dim(cube::YAXArray, candidates::Tuple)
    dim_names = Tuple(name.(cube.axes))
    for candidate in candidates
        if candidate in dim_names
            return candidate
        end
    end
    return nothing
end

function _align_longitude_reference(lon::AbstractVector, poly::AbstractMatrix)
    lon_aligned = Float64.(lon)
    poly_aligned = copy(poly)

    valid = .!isnan.(poly_aligned[1, :])
    if !any(valid)
        return lon_aligned, poly_aligned
    end

    if Base.maximum(lon_aligned) > 180 && any(poly_aligned[1, valid] .< 0)
        lon_aligned = map(x -> x >= 180 ? x - 360 : x, lon_aligned)
    elseif Base.maximum(lon_aligned) <= 180 && any(poly_aligned[1, valid] .> 180)
        poly_aligned[1, valid] = ifelse.(poly_aligned[1, valid] .> 180, poly_aligned[1, valid] .- 360, poly_aligned[1, valid])
    end

    return lon_aligned, poly_aligned
end

function spatialsubset(cube::YAXArray, poly)
    poly2 = _normalize_polygon(poly)

    lonsymbol = _find_spatial_dim(cube, (:longitude, :lon, :x, :rlon))
    latsymbol = _find_spatial_dim(cube, (:latitude, :lat, :y, :rlat))

    lonsymbol === nothing && error("Could not find longitude dimension in cube.")
    latsymbol === nothing && error("Could not find latitude dimension in cube.")

    lon = collect(lookup(cube, lonsymbol))
    lat = collect(lookup(cube, latsymbol))

    lon_mask, poly_mask = _align_longitude_reference(lon, poly2)

    longrid, latgrid = ndgrid(lon_mask, Float64.(lat))
    msk = inpolygrid(longrid, latgrid, poly_mask)

    I = findall(!isnan, msk)
    isempty(I) && error("Polygon does not overlap cube grid.")

    idlon = getindex.(I, 1)
    idlat = getindex.(I, 2)

    minXgrid = Base.minimum(idlon)
    maxXgrid = Base.maximum(idlon)
    minYgrid = Base.minimum(idlat)
    maxYgrid = Base.maximum(idlat)

    names = collect(name.(cube.axes))
    lonpos = findfirst(==(lonsymbol), names)
    latpos = findfirst(==(latsymbol), names)

    indices = Any[Colon() for _ in 1:ndims(cube)]
    indices[lonpos] = minXgrid:maxXgrid
    indices[latpos] = minYgrid:maxYgrid

    cubesub = view(cube, indices...)
    msksub = msk[minXgrid:maxXgrid, minYgrid:maxYgrid]

    maskshape = ntuple(i -> i == lonpos ? size(msksub, 1) : (i == latpos ? size(msksub, 2) : 1), ndims(cubesub))
    masked_data = cubesub.data .* reshape(Float32.(msksub), maskshape)

    return YAXArray(cubesub.axes, masked_data)
end

function spatialsubset(C::AbstractVector, poly)
    isempty(C) && error("Input container is empty.")
    C[1] isa YAXArray || error("Expected first element of input container to be a YAXArray.")

    out = copy(C)
    out[1] = spatialsubset(C[1], poly)
    return out
end

function spatialsubset(C::Tuple, poly)
    isempty(C) && error("Input container is empty.")
    C[1] isa YAXArray || error("Expected first element of input container to be a YAXArray.")

    return (spatialsubset(C[1], poly), Base.tail(C)...)
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

    minpoly = Base.minimum(poly[1, .!isnan.(poly[1, :])])
    maxpoly = Base.maximum(poly[1, .!isnan.(poly[1, :])])

    meridian = false
    if sign(minpoly)*sign(maxpoly) == sign(-1.0) # polygon crosses the meridian
       meridian = true
    end

    return meridian

end

function findmindist(p1::Tuple{Float64, Float64}, p2::Tuple{Array{Float64,1},Array{Float64,1}}; R=6372.8)

    dist = Array{Float64}(undef, length(p2[1]), 1)

    for idx = 1:length(p2[1])
        p3 = ( p2[1][idx], p2[2][idx] )
        dist[idx] = haversine(p1, p3, R)
    end

    return findmin(dist)

end


