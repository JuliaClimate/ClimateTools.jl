"""
    windnr(p, poly::Matrix)

Determines the winding number of a point and a polygon, i.e. how many
times a polygon winds around the point.

It follows Dan Sunday: http://geomalgorithms.com/a03-_inclusion.html.
"""
function windnr(p, poly::Matrix)
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
          (a closed poly)

Returns true if point has an odd winding number.  This should label
points as exterior which are inside outcrops.  See test for a test.

Author: Github "Mauro3" / "Mauro"
"""
function inpoly(p, poly::Matrix)
  return isodd(windnr(p, poly))
end
# inpoly(p, poly::Matrix) = isodd(windnr(p,poly))

"""
    [X, Y] = meshgrid(XV, YV)

This function creates a 2-D mesh-grid in a format consistent with Matlab's function meshgrid(). XV and YV are vectors.
"""
function meshgrid(XV, YV)

	# Number of columns ("width")
	num_cols = length(XV);

	# Number of rows ("height")
	num_rows = length(YV);

	# Replace the columns vector
	X = repmat(XV', num_rows, 1);

	# Replate the rows vector
	Y = repmat(YV, 1, num_cols);

	# Return X and Y
	return X, Y

end

"""
    inpolyvec(lon, lat, poly::AbstractArray{N,2} where N)

Used to test a vector of points. Columns should be consistent with polygon.


"""

function inpolyvec(lon, lat, poly::AbstractArray{N,2} where N)

    OUT = fill(NaN, size(lon, 1), size(lat, 1)) # grid mask

    # Convert longitude to 0-360 degrees_east
    lon[lon .< 0] += 360

    # Find number of polygons (separated by NaN values)
    polyidx = findn(isnan.(poly[1, :])) #poly start index
    npoly = length(polyidx) # number of polygons

    for p = 1:npoly # loop over each polygon
        # Build poly n
        if p != npoly
            polyn = poly[:, polyidx[p] + 1:polyidx[p + 1] - 1]
        elseif p == npoly
            polyn = poly[:, polyidx[p] + 1:end]
        end
        # Convert to 0-360 longitude
        if sum(polyn[1, :] .< 0) == length(polyn[1, :])
            polyn[1,:] += 360
        end

        # get_limits of polyn
        minlon = minimum(polyn[1, :])
        maxlon = maximum(polyn[1, :])
        minlat = minimum(polyn[2, :])
        maxlat = maximum(polyn[2, :])

        lonidx = findn((lon .<= maxlon) .& (lon .>= minlon))
        latidx = findn((lat .<= maxlat) .& (lat .>= minlat))


        for j in lonidx # = 1:size(lon, 1)
            for i in latidx # = 1:size(lat, 1) # loop over pair of points to test

                # if inpoly([lon[j], lat[i]], polyn)
                #     OUT[j, i] = 1.0
                # end
                if OUT[j, i] != 1.0 && inpoly([lon[j], lat[i]], polyn)
                    OUT[j, i] = 1.0
                end
            end
        end
    end
    return OUT

end


"""
    C = interp_climgrid(A::ClimGrid, B::ClimGrid)

This function interpolate ClimGrid A onto lat-lon grid of ClimGrid B,
where A, B and C are ClimGrid

"""
# TODO Add method where interp_climgrid(A::ClimGrid, B::lat/lon_grid)
function interp_climgrid(A::ClimGrid, B::ClimGrid)
    # ---------------------------------------
    # Get lat-lon information from ClimGrid B
    londest = B[1][Axis{:lon}][:]
    latdest = B[1][Axis{:lat}][:]

    if B.lonunits == "degrees_east" && sum(sign.(londest) .== -1) == length(londest) #indicate all negative longitude -> remap to 0-360 degrees
        londest = londest + 360.
        # recreate B ClimGrid and change lonunits to "degrees_west"
        timedest = B[1][Axis{:time}][:]
        axisdata = AxisArray(B[1].data, Axis{:time}(timedest), Axis{:lon}(londest), Axis{:lat}(latdest))
        B = ClimateTools.ClimGrid(axisdata, model = B.model, experiment = B.experiment, run = B.run, filename = B.filename, dataunits = B.dataunits, latunits = B.latunits, lonunits = "degrees_west", variable = B.variable, typeofvar = B.variable, typeofcal = B.typeofcal)
    end

    # Get lat-lon information from ClimGrid A
    lonorig = A[1][Axis{:lon}][:]
    latorig = A[1][Axis{:lat}][:]

    # -----------------------------------------
    # Get initial data and time from ClimGrid A
    dataorig = A[1].data
    timeorig = A[1][Axis{:time}][:] # the function will need to loop over time

    # Create lat-lon range consistent with length of data
    latorigrange = linspace(latorig[1], latorig[end], length(latorig))
    lonorigrange = linspace(lonorig[1], lonorig[end], length(lonorig))

    # ---------------------
    # Allocate output Array
    OUT = zeros(Float64, (length(timeorig), size(B.data, 2), size(B.data, 3)))

    # ------------------------
    # Interpolation
    for t = 1:length(timeorig)

        itp = interpolate(dataorig[t, :, :], BSpline(Linear()), OnCell())
        sitp = scale(itp, lonorigrange, latorigrange)

        for ilon = 1:length(londest)
            for ilat = 1:length(latdest)
                OUT[t, ilon, ilat] = sitp[londest[ilon], latdest[ilat]]
            end
        end

    end

    # -----------------------
    # Construct AxisArrays and ClimGrid struct from array OUT
    dataOut = AxisArray(OUT, Axis{:time}(timeorig), Axis{:lon}(londest), Axis{:lat}(latdest))

    C = ClimateTools.ClimGrid(dataOut, model = A.model, experiment = A.experiment, run = A.run, filename = A.filename, dataunits = A.dataunits, latunits = B.latunits, lonunits = B.lonunits, variable = A.variable, typeofvar = A.variable, typeofcal = A.typeofcal)
end

"""
    C = interp_climgrid(A::ClimGrid, londest::AbstractArray{N, 1} where N, latdest::AbstractArray{N, 1} where N)

This function interpolate ClimGrid A onto lat-lon grid defined by londest and latdest vector.

"""

function interp_climgrid(A::ClimGrid, londest::AbstractArray{N, 1} where N, latdest::AbstractArray{N, 1} where N)

    # Convert longitude to 0-360 degrees
    londest[londest .< 0.] += 360.

    # Get lat-lon information from ClimGrid A
    lonorig = A[1][Axis{:lon}][:]
    latorig = A[1][Axis{:lat}][:]

    # -----------------------------------------
    # Get initial data and time from ClimGrid A
    dataorig = A[1].data
    timeorig = A[1][Axis{:time}][:] # the function will need to loop over time

    # Create lat-lon range consistent with length of data
    latorigrange = linspace(latorig[1], latorig[end], length(latorig))
    lonorigrange = linspace(lonorig[1], lonorig[end], length(lonorig))

    # ---------------------
    # Allocate output Array
    OUT = zeros(Float64, (length(timeorig), length(londest), length(latdest)))

    # ------------------------
    # Interpolation
    for t = 1:length(timeorig)

        itp = interpolate(dataorig[t, :, :], BSpline(Linear()), OnCell())
        sitp = scale(itp, lonorigrange, latorigrange)

        for ilon = 1:length(londest)
            for ilat = 1:length(latdest)
                OUT[t, ilon, ilat] = sitp[londest[ilon], latdest[ilat]]
            end
        end

    end

    # -----------------------
    # Construct AxisArrays and ClimGrid struct from array OUT
    dataOut = AxisArray(OUT, Axis{:time}(timeorig), Axis{:lon}(londest), Axis{:lat}(latdest))

    C = ClimateTools.ClimGrid(dataOut, model = A.model, experiment = A.experiment, run = A.run, filename = A.filename, dataunits = A.dataunits, latunits = "degrees_north", lonunits = "degrees_east", variable = A.variable, typeofvar = A.variable, typeofcal = A.typeofcal)
end

"""
    applymask(A::AbstractArray{N, n}, mask::AbstractArray{N, n})

This function applies a mask on the array A. Return an AbstractArray{N, n}.

"""

function applymask(A::AbstractArray{N,4} where N, mask::AbstractArray{N, 2} where N)
    for t = 1:size(A, 1) # time axis
        for lev = 1:size(A, 4) #level axis
            tmp = A[t, :, :, lev]
            tmp .*= mask
            A[t, :, :, lev] = tmp
        end
    end
    return A
end

function applymask(A::AbstractArray{N,3} where N, mask::AbstractArray{N, 2} where N)
    for t = 1:size(A, 1) # time axis
        tmp = A[t, :, :]
        tmp .*= mask
        A[t, :, :] = tmp
    end
    return A
end

function applymask(A::AbstractArray{N,2} where N, mask::AbstractArray{N, 2} where N)
    A .*= mask
    return A
end

function applymask(A::AbstractArray{N,1} where N, mask::AbstractArray{N, 1} where N)
    A .*= mask
    return A
end
