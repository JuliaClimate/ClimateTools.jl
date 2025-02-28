"""
    C = griddata(A::ClimGrid, B::ClimGrid; min=[], max=[])

Interpolate `ClimGrid` A onto the lon-lat grid of `ClimGrid` B,
where A and B are `ClimGrid`.

Min and max optional keyword are used to constraint the results of the interpolation. For example, interpolating bounded fields can lead to unrealilstic values, such as negative precipitation. In that case, one would use min=0.0 to convert negative precipitation to 0.0.

"""
function griddata(A::ClimGrid, B::ClimGrid; method="linear", min=[], max=[])

    # ---------------------------------------
    # Get lat-lon information from ClimGrid B
    londest, latdest = ClimateTools.getgrids(B)

    # Get lat-lon information from ClimGrid A
    lonorig, latorig = ClimateTools.getgrids(A)
    points = hcat(lonorig[:], latorig[:])

    # -----------------------------------------
    # Get initial data and time from ClimGrid A
    dataorig = A[1].data
    timeorig = get_timevec(A)#[1][Axis{:time}][:] # the function will need to loop over time

    # ---------------------
    # Allocate output Array
    OUT = zeros(Float64, (size(B.data, 1), size(B.data, 2), length(timeorig)))

    # ------------------------
    # Interpolation
    # interp!(OUT, timeorig, dataorig, lonorig, latorig, londest, latdest, A.variable, frac=frac)
    griddata!(OUT, timeorig, dataorig, points, londest, latdest, method=method, msk=B.msk)

    if !isempty(min)
        OUT[OUT .<= min] .= min
    end

    if !isempty(max)
        OUT[OUT .>= max] .= max
    end

    # -----------------------
    # Construct AxisArrays and ClimGrid struct from array OUT
    latsymbol = Symbol(B.dimension_dict["lat"])
    lonsymbol = Symbol(B.dimension_dict["lon"])
    dataOut = AxisArray(OUT, Axis{lonsymbol}(B[1][Axis{lonsymbol}][:]), Axis{latsymbol}(B[1][Axis{latsymbol}][:]), Axis{:time}(timeorig))

    C = ClimateTools.ClimGrid(dataOut, longrid=B.longrid, latgrid=B.latgrid, msk=B.msk, grid_mapping=B.grid_mapping, dimension_dict=B.dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=A.experiment, run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits=B.latunits, lonunits=B.lonunits, variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)

end

"""
    C = griddata(A::ClimGrid, londest::AbstractArray{N, 1} where N, latdest::AbstractArray{N, 1} where N)

Interpolate `ClimGrid` A onto lat-lon grid defined by londest and latdest vector or array. If an array is provided, it is assumed that the grid is curvilinear (not a regular lon-lat grid) and the user needs to provide the dimension vector ("x" and "y") for such a grid.

"""
function griddata(A::ClimGrid, lon::AbstractArray{N, T} where N where T, lat::AbstractArray{N, T} where N where T; method="linear", dimx=[], dimy=[], min=[], max=[])

    # Get lat-lon information from ClimGrid A
    lonorig, latorig = getgrids(A)
    points = hcat(lonorig[:], latorig[:])

    # -----------------------------------------
    # Get initial data and time from ClimGrid A
    dataorig = A[1].data
    timeorig = A[1][Axis{:time}][:] # the function will need to loop over time

    if ndims(lon) == 1
        londest, latdest = ndgrid(lon, lat)
        dimx = lon
        dimy = lat
    elseif ndims(lon) == 2
        londest = lon
        latdest = lat
        @assert !isempty(dimx)
        @assert !isempty(dimy)
    else
        throw(error("Grid should be a vector or a grid"))
    end

    # ---------------------
    # Allocate output Array
    if ndims(lon) == 1
        OUT = zeros(Float64, (length(lon), length(lat), length(timeorig)))
    elseif ndims(lon) == 2
        OUT = zeros(Float64, (size(lon, 1), size(lon, 2), length(timeorig)))
    end

    # ------------------------
    # Interpolation
    griddata!(OUT, timeorig, dataorig, points, londest, latdest, method=method)

    if !isempty(min)
        OUT[OUT .<= min] .= min
    end

    if !isempty(max)
        OUT[OUT .>= max] .= max
    end

    # -----------------------
    # Construct AxisArrays and ClimGrid struct from array OUT
    dataOut = AxisArray(OUT, Axis{:x}(dimx), Axis{:y}(dimy), Axis{:time}(timeorig))
    msk = Array{Float64}(ones((size(OUT, 1), size(OUT, 2))))

    if ndims(lon) == 1
        grid_mapping = Dict(["grid_mapping_name" => "Regular_longitude_latitude"])
        dimension_dict = Dict(["lon" => "lon", "lat" => "lat"])
    elseif ndims(lon) == 2
        grid_mapping = Dict(["grid_mapping_name" => "Curvilinear_grid"])
        dimension_dict = Dict(["lon" => "x", "lat" => "y"])
    end

    C = ClimateTools.ClimGrid(dataOut, longrid=londest, latgrid=latdest, msk=msk, grid_mapping=grid_mapping, dimension_dict=dimension_dict, timeattrib=A.timeattrib, model=A.model, frequency=A.frequency, experiment=A.experiment, run=A.run, project=A.project, institute=A.institute, filename=A.filename, dataunits=A.dataunits, latunits="degrees_north", lonunits="degrees_east", variable=A.variable, typeofvar=A.typeofvar, typeofcal=A.typeofcal, varattribs=A.varattribs, globalattribs=A.globalattribs)

end

"""
    griddata!(OUT, timeorig, dataorig, points, londest, latdest, method, ;msk=[])
Interpolation of `dataorig` onto longitude grid `londest` and latitude grid `latdest`. Used internally by `regrid`.
"""
function griddata!(OUT, timeorig, dataorig, points, londest, latdest; method="linear", msk=[])

    p = Progress(length(timeorig), 5, "Regridding: ")
    for t = 1:length(timeorig)

        # Points values
        val = dataorig[:, :, t][:]

        # Call scipy griddata
        data_interp = scipy.griddata(points, val, (londest, latdest), method=method)

        # Apply mask from ClimGrid destination
        if !isempty(msk)
            OUT[:, :, t] .= data_interp .* msk
        else
            OUT[:, :, t] .= data_interp
        end

        next!(p)
    end
end

"""
    getgrids(C::ClimGrid)

Returns longitude and latitude grids of ClimGrid C.
"""
function getgrids(C::ClimGrid)
    longrid = C.longrid
    latgrid = C.latgrid
    return longrid, latgrid
end

"""
    getdims(C::ClimGrid)

Returns dimensions vectors of C
"""
function getdims(C::ClimGrid)
    latsymbol, lonsymbol = ClimateTools.getsymbols(C)
    x = C[1][Axis{lonsymbol}][:]
    y = C[1][Axis{latsymbol}][:]
    timevec = get_timevec(C)

    return x, y, timevec
end