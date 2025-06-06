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

# function applymask(C::ClimGrid, mask::AbstractArray{N, 1} where N)

#     A = C.data.data

#     modA = applymask(A, mask)

#     Anew = buildarrayinterface(modA, C)

#     return ClimGrid(Anew, longrid=C.longrid, latgrid=C.latgrid, msk=mask, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal="climatology", varattribs=C.varattribs, globalattribs=C.globalattribs)
# end

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

# """
#     ensemble_mean(C::ClimGrid...)

# Returns the Ensemble mean of ClimGrids C..
# """
# function ensemble_mean(C; skipnan=true)

#     # Create list of AxisArrays contained inside the ClimGrids
#     # Climref = C[1]
#     axisarrays = Array{Any}(undef, length(C))
#     # unitsClims = Array{Any}(undef, length(C))

#     for k = 1:length(C)
#         # datatmp[.!isnan.(datatmp)
#         axisarrays[k] = periodmean(C[k])[1]#[1][.!isnan.(C[k][1])], dims=3)
#         # unitsClims[k] = unit(C[k][1][1])
#     end

#     # if length(unique(unitsClims)) != 1
#     #     throw(error("ClimGrids needs to have the same physical units."))
#     # end

#     # ENSEMBLE MEAN
#     n = length(axisarrays) # number of members
#     dataout = sum(axisarrays) / n # ensemble mean

#     # # reapply unit
#     # dataout = [dataout][1]unitsClims[1]

#     # Build AxisArray
#     data_axis = ClimateTools.buildarrayinterface(dataout, C[1])


#     # Ensemble metadata
#     globalattribs = Dict()
#     globalattribs["history"] = "Ensemble mean"
#     globalattribs["models"] = ""
#     for k = 1:length(C)
#         globalattribs["models"] = string(globalattribs["models"], ", ", C[k].model)
#     end

#     # Return ClimGrid
#     return ClimGrid(data_axis, longrid=C[1].longrid, latgrid=C[1].latgrid, msk=C[1].msk, grid_mapping=C[1].grid_mapping, dimension_dict=C[1].dimension_dict, timeattrib=C[1].timeattrib, model=globalattribs["models"], frequency="Climatology ensemble mean", experiment="Multi-models ensemble", run="Multi-models ensemble", project="Multi-models ensemble", institute="Multi-models ensemble", filename="muliple_files", dataunits=C[1].dataunits, latunits=C[1].latunits, lonunits=C[1].lonunits, variable=C[1].variable, typeofvar=C[1].typeofvar, typeofcal="Climatology", varattribs=C[1].varattribs, globalattribs=globalattribs)

# end

# """
#     ensemble_std(C::ClimGrid...)

# Returns the Ensemble standard deviation of climatological means of ClimGrids C..
# """
# function ensemble_std(C; skipnan=true)

#     # Create list of AxisArrays contained inside the ClimGrids
#     # Climref = C[1]
#     axisarrays = Array{Any}(undef, length(C))
#     # unitsClims = Array{Any}(undef, length(C))

#     for k = 1:length(C)
#         # datatmp[.!isnan.(datatmp)
#         axisarrays[k] = periodmean(C[k])[1]#[1][.!isnan.(C[k][1])], dims=3)
#         # unitsClims[k] = unit(C[k][1][1])
#     end

#     # if typeof(unitsClims[1]) != Unitful.FreeUnits{(),NoDims,nothing}
#         # if length(unique(unitsClims)) != 1
#         #     throw(error("ClimGrids needs to have the same physical units."))
#         # end
#     # end

#     # ENSEMBLE STD
#     dataout = std(axisarrays)

#     # # reapply unit
#     # dataout = [dataout][1]unitsClims[1]

#     # Build AxisArray
#     data_axis = ClimateTools.buildarrayinterface(dataout, C[1])

#     # Ensemble metadata
#     globalattribs = Dict()
#     globalattribs["history"] = "Ensemble mean"
#     globalattribs["models"] = ""
#     for k = 1:length(C)
#         globalattribs["models"] = string(globalattribs["models"], ", ", C[k].model)
#     end

#     # Return ClimGrid
#     return ClimGrid(data_axis, longrid=C[1].longrid, latgrid=C[1].latgrid, msk=C[1].msk, grid_mapping=C[1].grid_mapping, dimension_dict=C[1].dimension_dict, timeattrib=C[1].timeattrib, model=globalattribs["models"], frequency="Climatology ensemble standard deviation", experiment="Multi-models ensemble", run="Multi-models ensemble", project="Multi-models ensemble", institute="Multi-models ensemble", filename="muliple_files", dataunits=C[1].dataunits, latunits=C[1].latunits, lonunits=C[1].lonunits, variable=C[1].variable, typeofvar=C[1].typeofvar, typeofcal="Climatology", varattribs=C[1].varattribs, globalattribs=globalattribs)

# end

# """
#     ensemble_max(C::ClimGrid...)

# Returns the Ensemble maximum of climatological means of ClimGrids C..
# """
# function ensemble_max(C; skipnan=true)

#     # Create list of AxisArrays contained inside the ClimGrids
#     # Climref = C[1]
#     axisarrays = Array{Any}(undef, length(C))
#     # unitsClims = Array{Any}(undef, length(C))

#     for k = 1:length(C)
#         # datatmp[.!isnan.(datatmp)
#         axisarrays[k] = periodmean(C[k])[1]#[1][.!isnan.(C[k][1])], dims=3)
#         # unitsClims[k] = unit(C[k][1][1])
#     end

#     # if typeof(unitsClims[1]) != Unitful.FreeUnits{(),NoDims,nothing}
#         # if length(unique(unitsClims)) != 1
#         #     throw(error("ClimGrids needs to have the same physical units."))
#         # end
#     # end

#     # ENSEMBLE MAX
#     dataout = maximum(axisarrays)

#     # # reapply unit
#     # dataout = [dataout][1]unitsClims[1]

#     # Build AxisArray
#     data_axis = ClimateTools.buildarrayinterface(dataout, C[1])

#     # Ensemble metadata
#     globalattribs = Dict()
#     globalattribs["history"] = "Ensemble mean"
#     globalattribs["models"] = ""
#     for k = 1:length(C)
#         globalattribs["models"] = string(globalattribs["models"], ", ", C[k].model)
#     end

#     # Return ClimGrid
#     return ClimGrid(data_axis, longrid=C[1].longrid, latgrid=C[1].latgrid, msk=C[1].msk, grid_mapping=C[1].grid_mapping, dimension_dict=C[1].dimension_dict, timeattrib=C[1].timeattrib, model=globalattribs["models"], frequency="Climatology ensemble maximum", experiment="Multi-models ensemble", run="Multi-models ensemble", project="Multi-models ensemble", institute="Multi-models ensemble", filename="muliple_files", dataunits=C[1].dataunits, latunits=C[1].latunits, lonunits=C[1].lonunits, variable=C[1].variable, typeofvar=C[1].typeofvar, typeofcal="Climatology", varattribs=C[1].varattribs, globalattribs=globalattribs)

# end

# """
#     ensemble_min(C::ClimGrid...)

# Returns the Ensemble minimum of climatological means of ClimGrids C..
# """
# function ensemble_min(C; skipnan=true)

#     # Create list of AxisArrays contained inside the ClimGrids
#     # Climref = C[1]
#     axisarrays = Array{Any}(undef, length(C))
#     # unitsClims = Array{Any}(undef, length(C))

#     for k = 1:length(C)
#         # datatmp[.!isnan.(datatmp)
#         axisarrays[k] = periodmean(C[k])[1]#[1][.!isnan.(C[k][1])], dims=3)
#         # unitsClims[k] = unit(C[k][1][1])
#     end

#     # if typeof(unitsClims[1]) != Unitful.FreeUnits{(),NoDims,nothing}
#         # if length(unique(unitsClims)) != 1
#         #     throw(error("ClimGrids needs to have the same physical units."))
#         # end
#     # end

#     # ENSEMBLE MIN
#     dataout = minimum(axisarrays)

#     # # reapply unit
#     # dataout = [dataout][1]unitsClims[1]

#     # Build AxisArray
#     data_axis = ClimateTools.buildarrayinterface(dataout, C[1])

#     # Ensemble metadata
#     globalattribs = Dict()
#     globalattribs["history"] = "Ensemble mean"
#     globalattribs["models"] = ""
#     for k = 1:length(C)
#         globalattribs["models"] = string(globalattribs["models"], ", ", C[k].model)
#     end

#     # Return ClimGrid
#     return ClimGrid(data_axis, longrid=C[1].longrid, latgrid=C[1].latgrid, msk=C[1].msk, grid_mapping=C[1].grid_mapping, dimension_dict=C[1].dimension_dict, timeattrib=C[1].timeattrib, model=globalattribs["models"], frequency="Climatology ensemble minimum", experiment="Multi-models ensemble", run="Multi-models ensemble", project="Multi-models ensemble", institute="Multi-models ensemble", filename="muliple_files", dataunits=C[1].dataunits, latunits=C[1].latunits, lonunits=C[1].lonunits, variable=C[1].variable, typeofvar=C[1].typeofvar, typeofcal="Climatology", varattribs=C[1].varattribs, globalattribs=globalattribs)

# end

function maximum(arrays::Array{Any})

    dataout = Array{Any}(undef, size(arrays[1], 1), size(arrays[1], 2))
    nbsims = length(arrays)

    I = size(arrays[1], 1)
    J = size(arrays[1], 2)

    for j = 1:J
        for i = 1:I
            datatmp = Array{Float64}(undef, nbsims)
            for k = 1:nbsims

                datatmp[k] = arrays[k][i, j]


            end
            dataout[i, j] = maximum(datatmp)
        end

    end
    return dataout

end


function minimum(arrays::Array{Any})

    dataout = Array{Any}(undef, size(arrays[1], 1), size(arrays[1], 2))
    nbsims = length(arrays)

    I = size(arrays[1], 1)
    J = size(arrays[1], 2)

    for j = 1:J
        for i = 1:I
            datatmp = Array{Float64}(undef, nbsims)
            for k = 1:nbsims

                datatmp[k] = arrays[k][i, j]


            end
            dataout[i, j] = minimum(datatmp)
        end

    end
    return dataout

end


# """
#     polyfit(C::ClimGrid)

# Returns an array of the polynomials functions of each grid points contained in ClimGrid C.
# """
# function polyfit(C::ClimGrid)

#     x = collect(1:length(C[1][Axis{:time}][:]))
#     dataout = Array{Polynomial{typeof(C[1][1,1,1])}}(undef, size(C[1], 1),size(C[1], 2))

#     # Reshaping for multi-threads
#     datain_rshp = reshape(C[1].data, size(C[1].data,1)*size(C[1].data,2), size(C[1].data,3))
#     dataout_rshp = reshape(dataout, size(dataout,1)*size(dataout,2),size(dataout,3))

#     # Threads.@threads for k = 1:size(datain_rshp,1)
#     for k = 1:size(datain_rshp,1)
#         y = datain_rshp[k,:]
#         polynomial = Polynomials.fit(x, y, 4)
#         polynomial[0] = 0.0
#         dataout_rshp[k] = polynomial

#     end

#     return dataout
# end

# """
#     polyval(C::ClimGrid, polynomial::Array{Poly{Float64}})

#     Returns a ClimGrid containing the values, as estimated from polynomial function polyn.
# """
# function polyval(C::ClimGrid, polynomial::Array{Polynomial{N},2} where N)

#     datain = C[1].data
#     dataout = Array{typeof(C[1][1,1,1])}(undef, (size(C[1], 1), size(C[1],2), size(C[1], 3)))

#     # Reshape
#     datain_rshp = reshape(datain, size(datain,1)*size(datain,2),size(datain,3))
#     dataout_rshp = reshape(dataout, size(dataout,1)*size(dataout,2),size(dataout,3))
#     polynomial_rshp = reshape(polynomial, size(polynomial,1)*size(polynomial,2))

#     # Threads.@threads for k = 1:size(datain_rshp,1)
#     for k = 1:size(datain_rshp,1)
#         val = polynomial_rshp[k].(datain_rshp[k,:])
#         dataout_rshp[k,:] = val
#     end

#     dataout2 = buildarrayinterface(dataout, C)

#     return ClimGrid(dataout2; longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
# end

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

# """
#     get_units(C::ClimGrid)
#
# Returns the adequate units of Unitful of ClimGrid C.
# """
# function get_units(C)
#     return unit(C[1][1,1,1])
# end

# """
#     uconvert(a::Units, C::ClimGrid)
#
# Convert a ClimGrid [`ClimateTools.ClimGrid`](@ref) to different units. The conversion will fail if the target units `a` have a different dimension than the dimension of the quantity `x`.
# """
#
# function uconvert(a::Units, C::ClimGrid)
#
#     dataconv = uconvert.(a, C[1])
#
#     dataaxis = buildarrayinterface(dataconv, C)
#
#     ClimGrid(dataaxis, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=string(unit(dataconv[1])), latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)
#
# end

"""
    replace_missing!(A::AbstractArray)

Replace missing value in array A with NaNs.
"""
function replace_missing!(A::AbstractArray)

    idx = findall(ismissing, A)
    if !isempty(idx)
        A[idx] .= NaN
    end
    return A
end

"""
    convert!(A::AbstractArray)

Convert the array A to a single type, removing the missing union.
"""
function convert!(A::AbstractArray, T)

    # Convert
    # T = typeof(A[1])
    B = Array{T}(A)
    return B

end

"""
    get_threshold(obsvec, refvec; thres=0.95)
"""
function get_threshold(obsvec, refvec; thres=0.95)
    return mean([quantile(obsvec[obsvec .>= 1.0],thres) quantile(refvec[refvec .>= 1.0], thres)])
end

# function timestep(C::ClimGrid, ts::Tuple)
#     D = temporalsubset(C, ts, ts)
# end

# function timestep(C::ClimGrid, ts::Int)
#     timevec = get_timevec(C)
#     d = timevec[ts]
#     t = (year(d), month(d), day(d))
#     D = temporalsubset(C, t, t)
# end

# function timestep(C::ClimGrid, d)
#     t = (year(d), month(d), day(d), hour(d), minute(d), second(d))
#     D = temporalsubset(C, t, t)
# end

# function get_max_clusters(x::Vector{Cluster})

    
#     dataout = Array{typeof(x[1].value[1])}(undef,0)

#     for i = eachindex(x)        
#         append!(dataout, x[i].value)
#     end


#     return dataout

# end

# function get_position_clusters(x::Vector{Cluster})

#     dataout = Array{Int}(undef, 0)

#     for i = eachindex(x)
#         append!(dataout, x[i].position)
#         #dataout[i] = x[i].position[1]
#     end


#     return dataout

# end


"""
    finitemean(x,y)

Calculate mean while omitting NaN, Inf, etc.
"""
finitemean(x) = mean(filter(x -> !isnan(x)&isfinite(x),x))
finitemean(x,y) = mapslices(finitemean,x,dims=y)

# """
#     periodmean(C::ClimGrid; startdate::Tuple, enddate::Tuple)

# Mean of array data over a given period.
# """
# function periodmean(C::ClimGrid; level=1, start_date::Tuple=(Inf, ), end_date::Tuple=(Inf,))

#     if start_date != (Inf,) || end_date != (Inf,)
#         C = temporalsubset(C, start_date, end_date)
#     end

#     datain = C.data.data

#     # Mean and squeeze
#     if ndims(datain) == 2
#         dataout = datain
#     elseif ndims(datain) == 3
#         if size(datain, 3) == 1 # already an average on single value
#             dataout = dropdims(datain, dims=3)
#         else
#             dataout = dropdims(finitemean(datain, 3), dims=3)
#         end
#     elseif ndims(datain) == 4
#         datain_lev = datain[:,:,level,:] # extract level
#         if size(datain_lev, 3) == 1
#             dataout = dropdims(datain_lev, dims=3)
#         else
#             dataout = dropdims(finitemean(datain_lev, 3), dims=3)
#         end
#     end

#     # Build output AxisArray
#     FD = buildarray_climato(C, dataout)

#     # Return climGrid type containing the indice
#     return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable="periodmean", typeofvar=C.typeofvar, typeofcal="climatology", varattribs=C.varattribs, globalattribs=C.globalattribs)
# end

# """
#     verticalmean(C::ClimGrid; startdate::Tuple, enddate::Tuple)

# Mean of array data over all vertical levels.
# """
# function verticalmean(C::ClimGrid)

#     datain = C.data.data

#     if ndims(datain) < 4
#         error("There is no vertical levels in the dataset")
#     end

#     if size(datain, 3) == 1 # Only one vertical level
#         dataout = dropdims(datain[:, :, :, :], dims = 3)
#     else
#         dataout = dropdims(finitemean(datain, 3), dims=3)
#     end

#     # Build output AxisArray
#     FD = buildarray_verticalmean(C, dataout)

#     # Return climGrid type containing the indice
#     return ClimGrid(FD, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable="periodmean", typeofvar=C.typeofvar, typeofcal="climatology", varattribs=C.varattribs, globalattribs=C.globalattribs)
# end

# function buildarray_climato(C::ClimGrid, dataout)
#     lonsymbol = Symbol(C.dimension_dict["lon"])
#     latsymbol = Symbol(C.dimension_dict["lat"])
#     if ndims(dataout) == 2 # original was a 3D field
#         FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}].val), Axis{latsymbol}(C[1][Axis{latsymbol}].val))
#     elseif ndims(dataout) == 3 # original was a 4D field
#         levsymbol = Symbol(C.dimension_dict["height"])
#         FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}].val), Axis{latsymbol}(C[1][Axis{latsymbol}].val), Axis{levsymbol}(C[1][Axis{levsymbol}].val))
#     end
#     return FD
# end

# function buildarray_verticalmean(C::ClimGrid, dataout)
#     lonsymbol = Symbol(C.dimension_dict["lon"])
#     latsymbol = Symbol(C.dimension_dict["lat"])
#     timesymbol = Symbol(C.dimension_dict["time"])

#     FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}].val), Axis{latsymbol}(C[1][Axis{latsymbol}].val), Axis{timesymbol}(C[1][Axis{timesymbol}].val))

#     return FD
# end

# function buildarrayinterface(axisArraytmp, A)
#     latsymbol = Symbol(A.dimension_dict["lat"])
#     lonsymbol = Symbol(A.dimension_dict["lon"])
#     if ndims(axisArraytmp) == 2
#         axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}].val), Axis{latsymbol}(A[1][Axis{latsymbol}].val))
#     elseif ndims(axisArraytmp) == 3
#         axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}].val), Axis{latsymbol}(A[1][Axis{latsymbol}].val), Axis{:time}(A[1][Axis{:time}].val))
#     elseif ndims(axisArraytmp) == 4
#         axisArray = AxisArray(axisArraytmp, Axis{lonsymbol}(A[1][Axis{lonsymbol}].val), Axis{latsymbol}(A[1][Axis{latsymbol}].val), Axis{:plev}(A[1][Axis{:plev}].val), Axis{:time}(A[1][Axis{:time}].val))
#     end
#     return axisArray
# end

# function buildarray_annual(C::ClimGrid, dataout, numYears)
#     lonsymbol = Symbol(C.dimension_dict["lon"])
#     latsymbol = Symbol(C.dimension_dict["lat"])
#     FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}].val), Axis{latsymbol}(C[1][Axis{latsymbol}].val), Axis{:time}(Dates.year.(DateTime.(numYears))))
#     return FD
# end

# function buildarray_resample(C::ClimGrid, dataout, newtime)
#     lonsymbol = Symbol(C.dimension_dict["lon"])
#     latsymbol = Symbol(C.dimension_dict["lat"])
#     FD = AxisArray(dataout, Axis{lonsymbol}(C[1][Axis{lonsymbol}].val), Axis{latsymbol}(C[1][Axis{latsymbol}].val), Axis{:time}(newtime))
#     return FD
# end


# """
#     function temporalsubset(C::ClimGrid, startdate::Date, enddate::Date)

# Returns the temporal subset of ClimGrid C. The temporal subset is defined by a start and end date.

# """
# function temporalsubset(C::ClimGrid, datebeg::Tuple, dateend::Tuple)

#     T = typeof(get_timevec(C)[1])
#     timeV = get_timevec(C)
#     idxtimebeg, idxtimeend = timeindex(timeV, datebeg, dateend, T)

#     # startdate = buildtimetype(datebeg, T)
#     # enddate = buildtimetype(dateend, T)

#     # some checkups
#     @argcheck idxtimebeg <= idxtimeend

#     dataOut = C[1][Axis{:time}(idxtimebeg:idxtimeend)]

#     # The following control ensure that a 1-timestep temporal subset returns a 3D Array with time information on the timestep. i.e. startdate == enddate
#     if ndims(dataOut) == 2
#         timeV = startdate
#         latsymbol = Symbol(C.dimension_dict["lat"])
#         lonsymbol = Symbol(C.dimension_dict["lon"])
#         data2 = fill(NaN, (size(dataOut,1), size(dataOut, 2), 1))
#         data2[:,:,1] = dataOut
#         dataOut = AxisArray(data2, Axis{lonsymbol}(C[1][Axis{lonsymbol}].val), Axis{latsymbol}(C[1][Axis{latsymbol}].val), Axis{:time}(C[1][Axis{:time}].val))

#     end

#     return ClimGrid(dataOut, longrid=C.longrid, latgrid=C.latgrid, msk=C.msk, grid_mapping=C.grid_mapping, dimension_dict=C.dimension_dict, timeattrib=C.timeattrib, model=C.model, frequency=C.frequency, experiment=C.experiment, run=C.run, project=C.project, institute=C.institute, filename=C.filename, dataunits=C.dataunits, latunits=C.latunits, lonunits=C.lonunits, variable=C.variable, typeofvar=C.typeofvar, typeofcal=C.typeofcal, varattribs=C.varattribs, globalattribs=C.globalattribs)

# end

# """
#     get_timevec(C::ClimGrid)

# Returns time vector of ClimGrid C.
# """
# get_timevec(C::ClimGrid) = C[1][Axis{:time}][:]

# """
#     buildtimetype(datetuple, f)

# Returns the adequate DateTime for temporal subsetting using DateType *f*
# """
# function buildtimetype(date_tuple, f)

#     if length(date_tuple) == 1
#         dateout = f(date_tuple[1], 01, 01)
#     elseif length(date_tuple) == 2
#         dateout = f(date_tuple[1], date_tuple[2], 01)
#     elseif length(date_tuple) == 3
#         dateout = f(date_tuple[1], date_tuple[2], date_tuple[3])
#     elseif length(date_tuple) == 4
#         dateout = f(date_tuple[1], date_tuple[2], date_tuple[3], date_tuple[4], 00, 00)
#     elseif length(date_tuple) == 5
#         dateout = f(date_tuple[1], date_tuple[2], date_tuple[3], date_tuple[4], date_tuple[5], 00)
#     elseif length(date_tuple) == 6
#         dateout = f(date_tuple[1], date_tuple[2], date_tuple[3], date_tuple[4], date_tuple[5], date_tuple[6])
#     end

#     return dateout
# end

# """
#     timeindex(timeVec, start_date, end_date, T)

# Return the index of time vector specified by start_date and end_date. T is the DateTime type (see NCDatasets.jl documentation).
# """
# function timeindex(timeV, datebeg::Tuple, dateend::Tuple, T)

#     # Start Date
#     if !isinf(datebeg[1])
#         # Build DateTime type
#         start_date = buildtimetype(datebeg, T)
#         # @argcheck start_date >= timeV[1]
#         idxtimebeg = findfirst(timeV .>= start_date)[1]
#     else
#         idxtimebeg = 1
#     end
#     # End date
#     if !isinf(dateend[1])
#         # Build DateTime type
#         end_date = buildtimetype(dateend, T)
#         # @argcheck end_date <= timeV[end]
#         idxtimeend = findlast(timeV .<= end_date)[1]
#     else
#         idxtimeend = length(timeV)
#     end

#     if !isinf(datebeg[1]) && !isinf(dateend[1])
#         @argcheck start_date <= end_date
#     end
#     return idxtimebeg, idxtimeend
# end


# """
#     getunit(C::ClimGrid)
#
# Returns the physical unit of ClimGrid C
# """
# function getunit(C::ClimGrid)
#     # Get unit
#     units_all = unique(unit.(C[1]))
#     return un = units_all[.!ismissing.(units_all)]
# end
#
# """
#     getunit(A::AbstractArray)
#
# Returns the physical unit of AbstractArray A
# """
# function getunit(A::AbstractArray)
#     # Get unit
#     units_all = unique(unit.(A))
#     return un = units_all[.!ismissing.(units_all)]
# end

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
