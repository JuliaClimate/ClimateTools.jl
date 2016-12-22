"""
Determines the winding number of a point and a polygon, i.e. how many
times a polygon winds around the point.

It follows Dan Sunday: http://geomalgorithms.com/a03-_inclusion.html.
"""
function windnr(p, poly::Matrix)
    @assert length(p)==2
    @assert poly[:,1]==poly[:,end]
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
Determines if a point is inside a polygon.

- p -- point (x,y) or [x,y]
- poly -- polygon vertices [x1 x2 ... xn x1
                            y1 y2 ... yn y1]
          (a closed poly)

Returns true if point has an odd winding number.  This should label
points as exterior which are inside outcrops.  See test for a test.
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
This function creates a boxcar averager with a window length of 3

function boxcar3(A::AbstractArray)

"""
function boxcar3(A::AbstractArray)
    out = similar(A)
    R = CartesianRange(size(A))
    I1, Iend = first(R), last(R)
    for I in R
        n, s = 0, zero(eltype(out))
        for J in CartesianRange(max(I1, I-I1), min(Iend, I+I1))
            s += A[J]
            n += 1
        end
        out[I] = s/n
    end
    out
end
