

"""
    regrid_cube(cube::YAXArray, target_grid::YAXArray; lonname_source = :longitude, latname_source=:latitude, lonname_dest=:longitude, latname_dest=:latitude) -> YAXArray

Regrids a `YAXArray` cube to a target grid using a linear interpolation between nodes.

# Arguments
- `cube::YAXArray`: The input data cube to be regridded.
- `target_grid::YAXArray`: The target grid to which the data cube will be regridded.
- `lonname_source::Symbol`: The longitude dimension name of the source dataset (defaults to :longitude).
- `latname_source::Symbol`: The latitude dimension name of the source dataset (defaults to :latitude).
- `lonname_dest::Symbol`: The longitude dimension name of the destination dataset (defaults to :longitude).
- `latname_dest::Symbol`: The latitude dimension name of the destination dataset (defaults to :latitude).

# Returns
- `YAXArray`: A new `YAXArray` that has been regridded to the target grid.


"""
function regrid_cube(source::YAXArray, dest::YAXArray; lonname_source = :longitude, latname_source=:latitude, lonname_dest=:longitude, latname_dest=:latitude)
    
    # subset "hardcodé" pour éviter effet de bords
    source_spatial = source[Dim{lonname_source}(Base.minimum(lookup(dest, lonname_dest))..Base.maximum(lookup(dest, lonname_dest))), Dim{latname_source}(Base.minimum(lookup(dest, latname_dest))..Base.maximum(lookup(dest, latname_dest)))]
    
    dest_spatial = dest[Dim{lonname_dest}(Base.minimum(lookup(source_spatial, lonname_source))..Base.maximum(lookup(source_spatial, lonname_source))), Dim{latname_dest}(Base.minimum(lookup(source_spatial,latname_source))..Base.maximum(lookup(source_spatial,latname_source)))]
    
    # Grille source
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

    # Grille destination    
    xg2, yg2 = collect(dest_spatial.longitude), collect(dest_spatial.latitude)
    Interpolations.deduplicate_knots!([xg2,yg2], move_knots = true)

    if Base.diff(xg2)[1] < 0.0 # vector needs to be monotonically increasing
        # dest_spatial = reverse(dest_spatial, dims=Dim{:longitude})
        xg2 = reverse(xg2)
    end

    if Base.diff(yg2)[1] < 0.0 # vector needs to be monotonically increasing
        # dest_spatial = reverse(dest_spatial, dims=Dim{:latitude})
        yg2 = reverse(yg2)
    end

    coords = Iterators.product(xg2,yg2)
    
    # Dimensions des calculs distribués
    indims = InDims(string(lonname_source),string(latname_source))
    outdims = OutDims(Dim{lonname_source}(xg2), Dim{latname_source}(yg2))

    # build_itp_spatial(xg1, yg1)

    return mapCube(regrid_cube, source_spatial; xg1=xg1, yg1=yg1, coords=coords, indims=indims, outdims=outdims)
    
end 

"""
    regrid_cube(xout, xin; xg1, yg1, coords)

"""
function regrid_cube(xout, xin; xg1, yg1, coords)

    # @show size(xin)
    
    # Interpolation structure            
    itp = interpolate((xg1,yg1), xin, Gridded(Linear()))    
    etpf = extrapolate(itp, Interpolations.Flat())
    
    return xout .= (c->etpf(c...)).(coords)
    
end