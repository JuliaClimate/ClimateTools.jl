
"""
    regrid_cube(source::YAXArray, dest::YAXArray; lonname_source = :longitude,
                latname_source = :latitude, lonname_dest = :longitude,
                latname_dest = :latitude) -> YAXArray

Regrids a `YAXArray` cube to a target grid using linear interpolation between nodes.
This implementation handles regular lon/lat grids. For curvilinear rotated-pole
grids use `regrid_rotated_curvilinear_to_regular` or `regrid_curvilinear_to_regular`.
"""
function regrid_cube(source::YAXArray, dest::YAXArray; lonname_source = :longitude, latname_source=:latitude, lonname_dest=:longitude, latname_dest=:latitude)

    # subset to avoid edge effects
    source_spatial = source[Dim{lonname_source}(Base.minimum(lookup(dest, lonname_dest))..Base.maximum(lookup(dest, lonname_dest))), Dim{latname_source}(Base.minimum(lookup(dest, latname_dest))..Base.maximum(lookup(dest,latname_dest)))]

    dest_spatial = dest[Dim{lonname_dest}(Base.minimum(lookup(source_spatial, lonname_source))..Base.maximum(lookup(source_spatial, lonname_source))), Dim{latname_dest}(Base.minimum(lookup(source_spatial,latname_source))..Base.maximum(lookup(source_spatial,latname_source)))]

    # Source grid
    xg1, yg1 = collect(source_spatial.longitude), collect(source_spatial.latitude)
    Interpolations.deduplicate_knots!([xg1,yg1], move_knots = true)
    
    if Base.diff(xg1)[1] < 0.0 # vector needs to be monotonically increasing
        source_spatial = reverse(source_spatial, dims=Dim{:longitude})
        xg1 = reverse(xg1)
    end
    
    if Base.diff(yg1)[1] < 0.0 # vector needs to be monotonically increasing
        source_spatial = reverse(source_spatial, dims=Dim{:latitude})
        yg1 = reverse(yg1)
    end
    
    # Destination grid
    xg2, yg2 = collect(dest_spatial.longitude), collect(dest_spatial.latitude)
    Interpolations.deduplicate_knots!([xg2,yg2], move_knots = true)
    
    if Base.diff(xg2)[1] < 0.0
        xg2 = reverse(xg2)
    end
    
    if Base.diff(yg2)[1] < 0.0
        yg2 = reverse(yg2)
    end
    
    coords = Iterators.product(xg2,yg2)
    
    # Dimensions for distributed calculations
    indims = InDims(string(lonname_source),string(latname_source))
    outdims = OutDims(Dim{lonname_source}(xg2), Dim{latname_source}(yg2))
    
    return mapCube(regrid_cube, source_spatial; xg1=xg1, yg1=yg1, coords=coords, indims=indims, outdims=outdims)
    
end

"""
    regrid_cube(xout, xin; xg1, yg1, coords)

Low-level per-chunk interpolation used by `mapCube`.
"""
function regrid_cube(xout, xin; xg1, yg1, coords)



    # Interpolation structure
    itp = interpolate((xg1,yg1), xin, Gridded(Linear()))
    etpf = extrapolate(itp, Interpolations.Flat())

    return xout .= (c->etpf(c...)).(coords)
end


# Utilities for rotated-pole curvilinear -> regular lon/lat regridding

"""
    idw_griddata(pts, vals, londest2, latdest2; k=8, p=2)

Simple inverse-distance-weighted (IDW) scattered interpolator implemented in pure Julia.
`pts` is an N×2 matrix with columns (lon, lat), `vals` is length-N, and
`londest2`, `latdest2` are destination 2D arrays. Returns a 2D Array of interpolated values.
"""
function idw_griddata(pts::AbstractMatrix, vals::AbstractVector, londest2, latdest2; k::Int=8, p::Real=2)
    nsrc = size(pts,1)
    nx, ny = size(londest2)
    out = Array{promote_type(eltype(vals),Float64)}(undef, nx, ny)

    # pre-extract source lon/lat vectors
    lon_src = pts[:,1]
    lat_src = pts[:,2]

    # Attempt to use NearestNeighbors.jl for k-NN queries if available (faster for large N)
    use_kdtree = false
    KDT = nothing
    try
        @eval begin
            import NearestNeighbors
        end
        use_kdtree = true
    catch
        use_kdtree = false
    end

    if use_kdtree
        # Build KDTree on lon/lat. Note: KDTree expects points as columns.
        # We'll adjust for longitude wrapping by duplicating longitudes shifted by ±360
        try
            # To handle longitude wrapping near the dateline, duplicate source points
            # shifted by ±360° in longitude. Keep track of original indices so we can
            # map KDTree results back to the original vals vector.
            lon_wrap = pts[:,1]
            lat_wrap = pts[:,2]
            pts_base = hcat(lon_wrap, lat_wrap)

            # create shifted copies
            pts_minus = hcat(lon_wrap .- 360.0, lat_wrap)
            pts_plus  = hcat(lon_wrap .+ 360.0, lat_wrap)

            # concatenated points (rows -> points); KDTree expects columns
            pts_all = vcat(pts_base, pts_minus, pts_plus)
            pts_mat = transpose(pts_all)

            # build KDTree
            kdt = NearestNeighbors.KDTree(pts_mat)

            # mapping from concatenated index back to original index
            nbase = size(pts_base,1)
            function map_index(idx)
                mod1(idx, nbase)
            end

            for j in 1:ny
                for i in 1:nx
                    xo = londest2[i,j]
                    yo = latdest2[i,j]

                    # query k nearest neighbors
                    ksel = min(k, nbase)
                    inds, dists = NearestNeighbors.knn(kdt, [xo, yo], ksel, true)
                    # dists are squared Euclidean distances; take sqrt
                    dsel = sqrt.(dists)

                    # exact match
                    if any(x->x < 1e-12, dsel)
                        rel = findfirst(x->x < 1e-12, dsel)
                        idx0 = map_index(inds[rel])
                        out[i,j] = vals[idx0]
                        continue
                    end

                    orig_inds = map(map_index, inds)
                    vsel = vals[orig_inds]
                    w = 1.0 ./ (dsel .^ p)
                    wsum = sum(w)
                    if wsum == 0.0
                        out[i,j] = NaN
                    else
                        out[i,j] = sum(w .* vsel) / wsum
                    end
                end
            end
        catch err
            @debug "KDTree IDW path failed, falling back to brute-force: $err"
            # fall back to brute-force below
            use_kdtree = false
        end
    end

    if !use_kdtree
        for j in 1:ny
            for i in 1:nx
                xo = londest2[i,j]
                yo = latdest2[i,j]

                # compute shortest longitude difference (wrap across ±180)
                dx = mod.(lon_src .- xo .+ 180.0, 360.0) .- 180.0
                dy = lat_src .- yo
                d = sqrt.(dx.^2 .+ dy.^2)

                # if any source point coincides exactly with the target, return its value
                idx0 = findfirst(x->iszero(x) || x < 1e-12, d)
                if idx0 !== nothing
                    out[i,j] = vals[idx0]
                    continue
                end

                # select k nearest neighbours (or all if fewer)
                ksel = min(k, nsrc)
                # partial sort to get indices of nearest k
                inds = partialsortperm(d, 1:ksel)
                dsel = d[inds]
                vsel = vals[inds]

                w = 1.0 ./ (dsel .^ p)
                wsum = sum(w)
                if wsum == 0.0
                    out[i,j] = NaN
                else
                    out[i,j] = sum(w .* vsel) / wsum
                end
            end
        end
    end

    return out
end

"""
    rotated_to_geographic(rlon, rlat, grid_north_longitude, grid_north_latitude)

Convert rotated-pole coordinates (rlon, rlat) to geographic longitude/latitude
using the CF rotated_pole convention. Inputs are in degrees. Returns (lon, lat)
arrays in degrees with the same shape as inputs.
"""
function rotated_to_geographic(rlon::AbstractArray, rlat::AbstractArray, grid_north_longitude::Real, grid_north_latitude::Real)
    # convert degrees -> radians
    deg2rad = pi/180.0
    φr = rlat .* deg2rad
    λr = rlon .* deg2rad
    β = grid_north_latitude * deg2rad
    λp = grid_north_longitude * deg2rad

    # CF rotated pole
    sinφ = sin(β) .* sin(φr) .+ cos(β) .* cos(φr) .* cos(λr)
    φ = asin.(sinφ)

    y = cos(φr) .* sin(λr)
    x = cos(β) .* sin(φr) .- sin(β) .* cos(φr) .* cos(λr)
    λ = λp .+ atan.(y, x)

    lat = φ ./ deg2rad
    lon = λ ./ deg2rad
    # normalize lon to [-180,180]
    lon = mod.(lon .+ 180.0, 360.0) .- 180.0

    return lon, lat
end


"""
    regrid_rotated_curvilinear_to_regular(lonrot, latrot, data, londest, latdest; grid_north_longitude, grid_north_latitude, method="linear")

Convert rotated-pole curvilinear (2D lon/lat) `lonrot`, `latrot` and source
`data` (2D) into geographic coordinates and interpolate onto the destination
grid (`londest`, `latdest`) using SciPy's griddata (via PyCall). `londest` and
`latdest` may be 1D vectors (regular grid) or 2D mesh arrays.
"""
function regrid_rotated_curvilinear_to_regular(lonrot::AbstractArray, latrot::AbstractArray, data::AbstractArray, londest, latdest; grid_north_longitude, grid_north_latitude, method::String="linear")
    @assert size(lonrot) == size(latrot) == size(data)

    # convert rotated -> geographic
    lont, latt = rotated_to_geographic(lonrot, latrot, grid_north_longitude, grid_north_latitude)

    # ensure destination grids are 2D arrays
    if ndims(londest) == 1 && ndims(latdest) == 1
        londest2, latdest2 = ndgrid(londest, latdest)
    else
        londest2 = londest
        latdest2 = latdest
    end

    # Try to detect if the source lon/lat are separable (rectilinear grid)
   (m, n) = size(lont)
    separable = false
    try
        # check if all columns of lont are equal (lon varies only with i)
        cols_equal = all([maximum(abs.(lont[:,j] .- lont[:,1])) < 1e-8 for j in 1:n])
        # check if all rows of latt are equal (lat varies only with j)
        rows_equal = all([maximum(abs.(latt[i,:] .- latt[1,:])) < 1e-8 for i in 1:m])
        if cols_equal && rows_equal
            separable = true
            xg1 = vec(lont[:,1])
            yg1 = vec(latt[1,:])
            # ensure increasing knots
            Interpolations.deduplicate_knots!([xg1, yg1], move_knots=true)
            if length(xg1) > 1 && Base.diff(xg1)[1] < 0.0
                xg1 = reverse(xg1)
                data = reverse(data, dims=1)
            end
            if length(yg1) > 1 && Base.diff(yg1)[1] < 0.0
                yg1 = reverse(yg1)
                data = reverse(data, dims=2)
            end

            # Gridded interpolation via Interpolations.jl
            itp = interpolate((xg1, yg1), data, Gridded(Linear()))
            etpf = extrapolate(itp, Interpolations.Flat())
            nx, ny = size(londest2)
            out = similar(londest2, Float64)
            for j in 1:ny
                for i in 1:nx
                    out[i,j] = etpf(londest2[i,j], latdest2[i,j])
                end
            end
            return out
        end
    catch err
        # any failure, fall back to IDW
        @debug "gridded interpolation detection failed: $err"
    end

    # Flatten points and values for IDW fallback
    pts = hcat(vec(lont), vec(latt))
    vals = vec(data)

    out = idw_griddata(pts, vals, londest2, latdest2; k=8, p=2)

    return out
end


"""
    regrid_curvilinear_to_regular(source, dest; grid_north_longitude, grid_north_latitude, method)

Wrapper that extracts rotated curvilinear lon/lat and data from common container
shapes used in this package (objects exposing `.longrid` / `.latgrid` and
`[1].data`) and calls `regrid_rotated_curvilinear_to_regular`.
"""
function regrid_curvilinear_to_regular(source, dest; grid_north_longitude, grid_north_latitude, method::String="linear")
    # try to extract rotated lon/lat and data from common container shapes
    lonrot = try
        source.longrid
    catch
        try
            source.longitude
        catch
            error("Could not find rotated longrid on source (expected .longrid or .longitude)")
        end
    end

    latrot = try
        source.latgrid
    catch
        try
            source.latitude
        catch
            error("Could not find rotated latgrid on source (expected .latgrid or .latitude)")
        end
    end

    data = try
        # ClimGrid-like: first AxisArray stored as source[1]
        source[1].data
    catch
        try
            Array(source)
        catch
            error("Could not extract data array from source")
        end
    end

    # destination grid
    londest = try
        dest.longrid
    catch
        try
            dest.longitude
        catch
            try
                dest[1][Dim{:longitude}][]
            catch
                error("Could not find destination grid on dest (expected .longrid/.longitude)")
            end
        end
    end
    latdest = try
        dest.latgrid
    catch
        try
            dest.latitude
        catch
            try
                dest[1][Dim{:latitude}][]
            catch
                error("Could not find destination grid on dest (expected .latgrid/.latitude)")
            end
        end
    end

    return regrid_rotated_curvilinear_to_regular(lonrot, latrot, data, londest, latdest; grid_north_longitude=grid_north_longitude, grid_north_latitude=grid_north_latitude, method=method)
end