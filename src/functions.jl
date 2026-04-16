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
    sz = length.(vs)
    return ntuple(n) do i
        shape = ntuple(j -> j == i ? sz[j] : 1, n)
        reps = ntuple(j -> j == i ? 1 : sz[j], n)
        repeat(reshape(vs[i], shape...), reps...)
    end
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
        minlon = Base.minimum(polyn[1, :])
        maxlon = Base.maximum(polyn[1, :])
        minlat = Base.minimum(polyn[2, :])
        maxlat = Base.maximum(polyn[2, :])

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

function inpolygrid(ds::Dataset, poly::AbstractArray{N,2} where N)

    londim, latdim = getdimnames(ds)

    longrid_flip = collect(ds["lon"].data)
    latgrid = collect(ds["lat"].data)


end 



"""
    applymask(A::AbstractArray{N, n}, mask::AbstractArray{N, n})

Applies a mask on the array A. Return an AbstractArray{N, n}.

"""
function applymask(A::AbstractArray{N,4} where N, mask::AbstractArray{N, 2} where N)

    T = typeof(A[.!ismissing.(A)][1])
    modA = Array{T}(undef, size(A))

    for t = 1:size(A, 4) # time axis
        for lev = 1:size(A, 3) #level axis
            tmp = A[:, :, lev, t]
            tmp .*= mask # TODO use multiple dispatch of applymask
            modA[:, :, lev, t] = tmp
        end
    end

    return modA
end

function applymask(A::AbstractArray{N,3} where N, mask::AbstractArray{N, 2} where N)

    T = typeof(A[.!ismissing.(A)][1])
    modA = Array{T}(undef, size(A))

    for t = 1:size(A, 3) # time axis
        tmp = A[:, :, t]
        tmp .*= mask # TODO use multiple dispatch of applymask
        modA[:, :, t] = tmp
    end

    return modA
end

function applymask(A::AbstractArray{N,2} where N, mask::AbstractArray{N, 2} where N)

    T = typeof(A[.!ismissing.(A)][1])
    modA = Array{T}(undef, size(A))

    @assert ndims(A) == ndims(mask)
    modA = A .* mask

    return modA
end

function applymask(A::AbstractArray{N,1} where N, mask::AbstractArray{N, 1} where N)

    T = typeof(A[.!ismissing.(A)][1])
    modA = Array{T}(undef, size(A))

    modA = A .* mask

    return modA
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


"""
    extension(url::String)

Returns the file extension of *url*.
"""

function extension(url::String)
    try
        return match(r"\.[A-Za-z0-9]+$", url).match
    catch
        return ""
    end
end



"""
    get_threshold(obsvec, refvec; thres=0.95)
"""
function get_threshold(obsvec, refvec; thres=0.95)
    return mean([quantile(obsvec[obsvec .>= 1.0],thres) quantile(refvec[refvec .>= 1.0], thres)])
end

"""
    finitemean(x,y)

Calculate mean while omitting NaN, Inf, etc.
"""
finitemean(x) = mean(filter(x -> !isnan(x)&isfinite(x),x))
finitemean(x,y) = mapslices(finitemean,x,dims=y)



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
