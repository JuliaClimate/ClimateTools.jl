"""
    windnr(p, poly::Matrix)

Determines the winding number of a point and a polygon, i.e. how many
times a polygon winds around the point.

It follows Dan Sunday: http://geomalgorithms.com/a03-_inclusion.html.
"""
function windnr(p, poly)
    @assert length(p) == 2
    @assert poly[:, 1] == poly[:, end]
    # Loop over edges
    px = p[1]
    py = p[2]
    wn = 0 # winding number, inside if >0
    len = length(poly)-2
    @inbounds for k = 1:2:len #size(poly,2)-1 # @fastmath makes it worse
        # ex1,ey1,ex2,ey2 = poly[k:k+3] # is slow...
        ex1, ey1, ex2, ey2 = poly[k], poly[k+1], poly[k+2], poly[k+3]
        # Reject edge if p totally above or below p:
        if (ey1>py && ey2>py) || (ey1<py && ey2<py)
            continue
        end
        # Check edges intersecting a horizontal ray:
        orient = leftorright(px,py,ex1,ey1,ex2,ey2)
        if up(ey1, ey2)
            if orient == -1; wn -= 1 end # p strictly left of e
        elseif down(ey1, ey2)
            if orient==1;  wn += 1 end # p strictly right of e
        end
    end
    return wn
end

up(ey1,ey2)   = ey1 < ey2
down(ey1,ey2) = ey1 > ey2

function leftorright(px,py,ex1,ey1,ex2,ey2)
    # returns:
    # -1 if on left of line
    #  0 if on line
    #  1 if on right of line
    vx,vy = ex2-ex1, ey2-ey1
    px,py =  px-ex1,  py-ey1
    sign(px * vy - py * vx)
end

"""
    inpoly(p, poly::Matrix)

Determines if a point is inside a polygon.

- p -- point (x,y) or [x,y]
- poly -- polygon vertices [x1 x2 ... xn x1
                            y1 y2 ... yn y1]


Returns true if point has an odd winding number.  This should label
points as exterior which are inside outcrops.  See test for a test.

Author: Github "Mauro3" / "Mauro"
"""
function inpoly(p, poly)
  return isodd(windnr(p, poly))
end
# inpoly(p, poly::Matrix) = isodd(windnr(p,poly))

"""
    X, Y = ndgrid(XV, YV)

This function creates a 2-D mesh-grid in a format consistent with Matlab's function ndgrid(). XV and YV are vectors.
"""
ndgrid(v::AbstractVector) = copy(v)

function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where T
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repeat(v1, 1, n), repeat(v2, m, 1))
end

function ndgrid(vs::AbstractVector{T}...) where T
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Array{T}(sz), n)
    s = 1
    for i=1:n
        a = out[i]::Array
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end

meshgrid(v::AbstractVector) = meshgrid(v, v)

"""
    X, Y = meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})

This function creates a 2-D mesh-grid in a format consistent with Matlab's function meshgrid(). XV and YV are vectors.
"""
function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where T
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repeat(vx, m, 1), repeat(vy, 1, n))
end

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T},
    vz::AbstractVector{T}) where T
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

"""
    inpolygrid(lon, lat, poly::AbstractArray{N,2} where N)

Used to test a grid of points. Returns a mask of ones and NaNs of the same size as lon and lat.

"""
function inpolygrid(lon::AbstractArray{N, 2} where N, lat::AbstractArray{N,2} where N, poly::AbstractArray{N,2} where N)

    OUT = fill(NaN, size(lon)) # grid mask

    # Find number of polygons (separated by NaN values)
    polyidx = Base.findall(isnan, poly[1,:])
    npoly = length(polyidx) # number of polygons

    for p = 1:npoly # loop over each polygon
        # Build poly n
        if p != npoly
            polyn = poly[:, polyidx[p] + 1:polyidx[p + 1] - 1]
        elseif p == npoly
            polyn = poly[:, polyidx[p] + 1:end]
        end

        # get_limits of polyn
        minlon = minimum(polyn[1, :])
        maxlon = maximum(polyn[1, :])
        minlat = minimum(polyn[2, :])
        maxlat = maximum(polyn[2, :])

        begin
            I = Base.findall((lon .<= maxlon) .& (lon .>= minlon) .& (lat .>= minlat) .& (lat .<= maxlat))
            idx, idy = (getindex.(I, 1), getindex.(I, 2))
        end

        for (ix, iy) in zip(idx, idy)
            # @show lon[ix, iy], lat[ix, iy]
            if OUT[ix, iy] != 1.0 && inpoly([lon[ix, iy], lat[ix, iy]], polyn)
                OUT[ix, iy] = 1.0
            end
        end

    end
    return OUT

end

# TODO define interpolation for 4D grid

"""
    C = regrid(A::ClimGrid, B::ClimGrid; method="linear", min=[], max=[])

Interpolate `ClimGrid` A onto the lon-lat grid of `ClimGrid` B,
where A and B are `ClimGrid`. Available methods for interpolation are "linear" (default), "nearest" and "cubic".

Min and max optional keyword are used to constraint the results of the interpolation. For example, interpolating bounded fields can lead to unrealilstic values, such as negative precipitation. In that case, one would use min=0.0 to convert negative precipitation to 0.0.

"""
function regrid(A::ClimGrid, B::ClimGrid; method::String="linear", min=[], max=[])

    # ---------------------------------------
    # Get lat-lon information from ClimGrid B
    londest, latdest = getgrids(B)

    # Get lat-lon information from ClimGrid A
    lonorig, latorig = getgrids(A)
    points = hcat(lonorig[:], latorig[:])

    # -----------------------------------------
    # Get initial data and time from ClimGrid A
    dataorig = view(A[1].data,:, :, :)
    timeorig = A[1][Axis{:time}][:] # the function will need to loop over time

    # ---------------------
    # Allocate output Array
    OUT = zeros(Float64, (size(B.data, 1), size(B.data, 2), length(timeorig)))

    # ------------------------
    # Interpolation
    interp!(OUT, timeorig, dataorig, points, londest, latdest, method, msk=B.msk)

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
    C = regrid(A::ClimGrid, londest::AbstractArray{N, 1} where N, latdest::AbstractArray{N, 1} where N)A

Interpolate `ClimGrid` A onto lat-lon grid defined by londest and latdest vector or array. If an array is provided, it is assumed that the grid is curvilinear (not a regular lon-lat grid) and the user needs to provide the dimension vector ("x" and "y") for such a grid.

"""
function regrid(A::ClimGrid, lon::AbstractArray{N, T} where N where T, lat::AbstractArray{N, T} where N where T; dimx=[], dimy=[], method::String="linear", min=[], max=[])

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
    interp!(OUT, timeorig, dataorig, points, londest, latdest, method)

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
    interp!(OUT, timeorig, dataorig, points, londest, latdest, method, ;msk=[])

Interpolation of `dataorig` onto longitude grid `londest` and latitude grid `latdest`. Used internally by `regrid`.
"""
function interp!(OUT, timeorig, dataorig, points, londest, latdest, method, ;msk=[])

    p = Progress(length(timeorig), 5, "Regridding: ")
    for t = 1:length(timeorig)

        # Points values
        val = dataorig[:, :, t][:]

        # Call scipy griddata
        data_interp = scipy[:griddata](points, val, (londest, latdest), method=method)

        # Apply mask from ClimGrid destination
        if !isempty(msk)
            OUT[:, :, t] = data_interp .* msk
        else
            OUT[:, :, t] = data_interp
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

"""
    get_timevec(C::ClimGrid)

Returns time vector of ClimGrid C.
"""
get_timevec(C::ClimGrid) = C[1][Axis{:time}][:]

"""
    applymask(A::AbstractArray{N, n}, mask::AbstractArray{N, n})

Applies a mask on the array A. Return an AbstractArray{N, n}.

"""
function applymask(A::AbstractArray{N,4} where N, mask::AbstractArray{N, 2} where N)
    for t = 1:size(A, 4) # time axis
        for lev = 1:size(A, 3) #level axis
            tmp = A[:, :, lev, t]
            tmp .*= mask # TODO use multiple dispatch of applymask
            A[:, :, lev, t] = tmp
        end
    end
    return A
end

function applymask(A::AbstractArray{N,3} where N, mask::AbstractArray{N, 2} where N)
    for t = 1:size(A, 3) # time axis
        tmp = A[:, :, t]
        tmp .*= mask # TODO use multiple dispatch of applymask
        A[:, :, t] = tmp
    end
    return A
end

function applymask(A::AbstractArray{N,2} where N, mask::AbstractArray{N, 2} where N)
    @assert ndims(A) == ndims(mask)
    A .*= mask
    return A
end

function applymask(A::AbstractArray{N,1} where N, mask::AbstractArray{N, 1} where N)
    A .*= mask
    return A
end

macro isdefined(var)
    quote
        try
            local _ = $(esc(var))
            true
        catch err
            isa(err, UndefVarError) ? false : rethrow(err)
        end
    end
end


function permute_west_east(data::AbstractArray{N,T} where N where T, longrid)#iwest, ieast)

    dataout = similar(data)
    ieast = longrid .>= 0.0
    iwest = longrid .< 0.0


    if ndims(data) == 2

        datatmp = data[:, :]
        # datawest = reshape(datatmp[iwest], :, size(datatmp, 2))
        # dataeast = reshape(datatmp[ieast], :, size(datatmp, 2))
        dataout[:, :] = permute_west_east2D(datatmp, iwest, ieast)


    elseif ndims(data) == 3

        for t = 1:size(data, 3)
            datatmp = data[:, :, t]
            # datawest = reshape(datatmp[iwest], :, size(datatmp, 2))
            # dataeast = reshape(datatmp[ieast], :, size(datatmp, 2))
            dataout[:, :, t] = permute_west_east2D(datatmp, iwest, ieast) #vcat(datawest, dataeast)
        end

    elseif ndims(data) == 4

        for t = 1:size(data, 4)
            for l = 1:size(data, 3)
                datatmp = data[:, :, l, t]
                # datawest = reshape(datatmp[iwest], :, size(datatmp, 2))
                # dataeast = reshape(datatmp[ieast], :, size(datatmp, 2))
                dataout[:, :, l, t] = permute_west_east2D(datatmp, iwest, ieast) #vcat(datawest, dataeast)
            end
        end

    end

    return dataout

end

function permute_west_east2D(data::AbstractArray{N,2} where N, iwest, ieast)

    datawest = reshape(data[iwest], :, size(data, 2))
    dataeast = reshape(data[ieast], :, size(data, 2))
    return vcat(datawest, dataeast)

end

function permute_east_west2D(data::AbstractArray{N,2} where N, iwest, ieast)

    datawest = reshape(data[iwest], :, size(data, 2))
    dataeast = reshape(data[ieast], :, size(data, 2))
    return vcat(dataeast, datawest)

end

"""
    ensemble_mean(C::ClimGrid...)

Returns the Ensemble mean of ClimGrids C..
"""
function ensemble_mean(C; skipnan=true)

    # Create list of AxisArrays contained inside the ClimGrids
    # Climref = C[1]
    axisarrays = Array{ClimGrid}(undef, length(C))

    for k = 1:length(C)
        # datatmp[.!isnan.(datatmp)
        axisarrays[k] = periodmean(C[k])#[1][.!isnan.(C[k][1])], dims=3)
    end

    # ENSEMBLE MEAN
    n = length(axisarrays) # number of members
    dataout = sum(axisarrays) / n # ensemble mean

    # Ensemble metadata
    globalattribs = Dict()
    globalattribs["history"] = "Ensemble mean"
    globalattribs["models"] = ""
    for k = 1:length(C)
        globalattribs["models"] = string(globalattribs["models"], ", ", C[k].model)
    end

    # Return ClimGrid
    return ClimGrid(dataout[1], longrid=C[1].longrid, latgrid=C[1].latgrid, msk=C[1].msk, grid_mapping=C[1].grid_mapping, dimension_dict=C[1].dimension_dict, timeattrib=C[1].timeattrib, model=globalattribs["models"], frequency="Climatology", experiment="Multi-models ensemble", run="Multi-models ensemble", project="Multi-models ensemble", institute="Multi-models ensemble", filename="muliple_files", dataunits=C[1].dataunits, latunits=C[1].latunits, lonunits=C[1].lonunits, variable=C[1].variable, typeofvar=C[1].typeofvar, typeofcal="Climatology", varattribs=C[1].varattribs, globalattribs=globalattribs)

end


"""
    polyfit(C::ClimGrid)

Returns an array of the polynomials functions of each grid points contained in ClimGrid C.
"""
function polyfit(C::ClimGrid)
    x = 1:length(C[1][Axis{:time}][:])
    dataout = Array{Polynomials.Poly{Float64}}(undef, size(C[1], 1),size(C[1], 2))

    # Reshaping for multi-threads
    datain_rshp = reshape(C[1].data, size(C[1].data,1)*size(C[1].data,2), size(C[1].data,3))
    dataout_rshp = reshape(dataout, size(dataout,1)*size(dataout,2),size(dataout,3))

    Threads.@threads for k = 1:size(datain_rshp,1)
        y = datain_rshp[k,:]
        polynomial = Polynomials.polyfit(x, y, 4)
        polynomial[0] = 0.0
        dataout_rshp[k] = polynomial

    end



    # for k = 1:size(C[1], 2)
    #     Threads.@threads for j = 1:size(C[1], 1)
    #         y = C[1][j , k, :].data
    #         polynomial = Polynomials.polyfit(x, y, 4)
    #         polynomial[0] = 0.0
    #         dataout[j,k] = polynomial
    #     end
    # end
    return dataout
end

"""
    polyval(C::ClimGrid, polynomial::Array{Poly{Float64}})

    Returns a ClimGrid containing the values, as estimated from polynomial function polyn.
"""
function polyval(C::ClimGrid, polynomial::Array{Poly{Float64},2})
    datain = C[1].data
    dataout = fill(NaN, (size(C[1], 1), size(C[1],2), size(C[1], 3)))::Array{N, T} where N where T

    # Reshape
    datain_rshp = reshape(datain, size(datain,1)*size(datain,2),size(datain,3))
    dataout_rshp = reshape(dataout, size(dataout,1)*size(dataout,2),size(dataout,3))
    polynomial_rshp = reshape(polynomial, size(polynomial,1)*size(polynomial,2))

    Threads.@threads for k = 1:size(datain_rshp,1)
        val = polynomial_rshp[k](datain_rshp[k,:])
        dataout_rshp[k,:] = val
    end

    #
    # for k = 1:size(C[1], 2)
    #     Threads.@threads for j = 1:size(C[1], 1)
    #         val = polynomial[j,k](datain[j,k,:])
    #         dataout[j,k,:] = val
    #     end
    # end

    dataout2 = buildarrayinterface(dataout, C)

    return ClimGrid(dataout2; longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
end

function extension(url::String)
    try
        return match(r"\.[A-Za-z0-9]+$", url).match
    catch
        return ""
    end
end

# function rot2lonlat(lon, lat, SP_lon, SP_lat; northpole = true)
#
#     # Copyright (c) 2013, Simon Funder
#     # All rights reserved.
#     #
#     # Redistribution and use in source and binary forms, with or without
#     # modification, are permitted provided that the following conditions are
#     # met:
#     #
#     #     * Redistributions of source code must retain the above copyright
#     #       notice, this list of conditions and the following disclaimer.
#     #     * Redistributions in binary form must reproduce the above copyright
#     #       notice, this list of conditions and the following disclaimer in
#     #       the documentation and/or other materials provided with the distribution
#     #
#     # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#     # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#     # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#     # ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
#     # LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#     # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#     # SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#     # INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#     # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#     # ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#     # POSSIBILITY OF SUCH DAMAGE.
#
#
#
#     # Convert degrees to radians
#     lon = (lon*π) ./ 180.0
#     lat = (lat*π) ./ 180.0
#
#     if northpole
#         SP_lon = SP_lon - 180.0
#         SP_lat = -SP_lat
#
#         θ = 90.0 + SP_lat # Rotation around y-axis
#         # θ = NP_lat # Rotation around y-axis
#         ϕ = SP_lon # Rotation around z-axis
#     else
#
#         θ = 90.0 + SP_lat # Rotation around y-axis
#         # θ = NP_lat # Rotation around y-axis
#         ϕ = SP_lon # Rotation around z-axis
#     end
#
#     ϕ = (ϕ * π) / 180.0 # Convert degrees to radians
#     θ = (θ * π) / 180.0
#
#     # Convert from spherical to cartesian coords
#     x = cos(lon) * cos(lat)
#     y = sin(lon) * cos(lat)
#     z = sin(lat)
#
#     ϕ = -ϕ
#     θ = -θ
#
#     x_new = cos(θ)*cos(ϕ)*x + sin(ϕ)*y + sin(θ).*cos(ϕ)*z
#     y_new = -cos(θ)*sin(ϕ)*x + cos(ϕ)*y - sin(θ)*sin(ϕ)*z
#     z_new = -sin(θ)*x + cos(θ)*z
#
#     # Convert cartesian back to spherical coordinates
#     lon_new = atan2(y_new, x_new)
#     lat_new = asin(z_new)
#
#     # Convert radians back to degrees
#     lon_new = (lon_new * 180.0) / π
#     lat_new = (lat_new * 180.0) / π
#
#     return lon_new, lat_new
#
#
# end

# # example for rot2lonlat
# SP_lon2 = 18
# SP_lat2 = -39.3
#
# grid_in = [[12; 12; 12] [55; 54; 53]]
# rot2lonlat(grid_in[:, 1], grid_in[:, 2], SP_lon2, SP_lat2, northpole=false)
