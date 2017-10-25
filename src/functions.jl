"""
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
This function creates a 2-D mesh-grid in a format consistent with Matlab's function meshgrid()

[X, Y] = meshgrid(XV, YV)

where XV and YV are vectors.
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
Used to test a vector of points. Columns should be consistent with polygon.
inpolyvec(pts::Array{Float64,2}, poly::Array{Float64,2})

"""

function inpolyvec(pts::Array{Float64,2}, poly::Array{Float64,2})

    OUT = fill(NaN, size(pts, 1))

    for i = 1:size(pts, 1) # loop over pair of points to test
        if inpoly([pts[i, 1], pts[i, 2]], poly)
            OUT[i] = 1.0
        end

    end
    return OUT

end

function inpolyvec(lon, lat, poly::Array{Float64,2})

    OUT = fill(NaN, size(lon, 1), size(lat, 1))

    for i = 1:size(lat, 1) # loop over pair of points to test
        for j = 1:size(lon, 1)
            if inpoly([lon[j], lat[i]], poly)
                OUT[j, i] = 1.0
            end
        end
    end
    return OUT

end


"""
This function interpolate ClimGrid A onto lat-lon grid of ClimGrid B

C = interpolate(A, B)

where A, B and C are ClimGrid

"""

function interp_climgrid(A::ClimGrid, B::ClimGrid)
    # ---------------------------------------
    # Get lat-lon information from ClimGrid B
    londest = B[1][Axis{:lon}][:]
    latdest = B[1][Axis{:lat}][:]

    if B.lonunits == "degrees_east" && sum(sign(londest) .== -1) == length(londest) #indicate all negative longitude -> remap to 0-360 degrees
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

    # Ti = scipy[:griddata]([lonorig, latorig], dataorig[1, :, :], [X, Y], method = "cubic")

    # ------------------------
    # Interpolation
    for t = 1:length(timeorig)
        # spline = Spline2D(lonorig, latorig, dataorig[t, :, :])
        # OUT[t, :, :] = evalgrid(spline, X[1, 500:550], Y[460:-1:400, 1])

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

function applymask(A::Array{Float64,4}, mask)
    for t = 1:size(A, 4)
        for lev = 1:size(A, 3)
            tmp = A[:, :, lev, t]
            tmp .*= mask
            A[:, :, lev, t] = tmp
        end
    end
    return A
end

function applymask(A::Array{Float64,3}, mask)
    for t = 1:size(A, 3)
        tmp = A[:, :, t]
        tmp .*= mask
        A[:, :, t] = tmp
    end
    return A
end

function applymask(A::Array{Float64,2}, mask)
    A .*= mask
    return A
end

function applymask(A::Array{Float32,4}, mask)
    for t = 1:size(A, 4)
        for lev = 1:size(A, 3)
            tmp = A[:, :, lev, t]
            tmp .*= mask
            A[:, :, lev, t] = tmp
        end
    end
    return A
end

function applymask(A::Array{Float32,3}, mask)
    for t = 1:size(A, 3)
        tmp = A[:, :, t]
        tmp .*= mask
        A[:, :, t] = tmp
    end
    return A
end

function applymask(A::Array{Float32,2}, mask)    
    A .*= mask
    return A
end
