function rlevels_cube(xout, xin; threshold=nothing, rlevels = [1, 2, 5, 10, 25, 50, 100, 1000], minimalvalue=1.0)

    if all(ismissing, xin)
        xout .= missing
        return    
    end

    if isnothing(threshold)
        dataforquantile = xin[xin .> minimalvalue]
        if  !isempty(dataforquantile)
            threshold = quantile(dataforquantile, 0.92)
            exceedances = xin[xin .> threshold] .- threshold
            model = Extremes.gpfit(exceedances)

            nobs = size(xin,1)        

            nobsperblock = 365

            r_h = returnlevel.(model, threshold, nobs, nobsperblock, rlevels)

            xout .= [x.value[1] for x in r_h]

        else       
            xout .= missing
            return
        end
        
    end

end


"""
    rlevels_cube(ds::YAXArray; rlevels = [1, 2, 5, 10, 25, 50, 100, 1000], threshold=nothing, minimalvalue=1.0, latsouth_pixel=0, latnorth_pixel=0, lonwest_pixel=0, loneast_pixel=0)

Compute the return levels for a dataset `ds` along the time dimension.

# Arguments
- `ds`: The input dataset.
- `rlevels`: The return levels to compute. Default is [1, 2, 5, 10, 25, 50, 100, 1000].
- `threshold`: The threshold to use for the Generalized Pareto fit. If not provided, the 0.92 quantile is used.
- `minimalvalue`: The minimal value to consider in the dataset. Default is 1.0.
- `latsouth_pixel`: The south pixel to consider. Default is 0.
- `latnorth_pixel`: The north pixel to consider. Default is 0.
- `lonwest_pixel`: The west pixel to consider. Default is 0.
- `loneast_pixel`: The east pixel to consider. Default is 0.

# Returns
The return levels of the input dataset `ds` along the time dimension.

"""
function rlevels_cube(ds::YAXArray; threshold=nothing, rlevels = [1, 2, 5, 10, 25, 50, 100, 1000], minimalvalue=1.0, latsouth_pixel=0, latnorth_pixel=0, lonwest_pixel=0, loneast_pixel=0, lonname="longitude", latname="latitude")

    # Dimensions
    # indims = InDims(MovingWindow(latname, latsouth_pixel,latnorth_pixel), MovingWindow(lonname, lonwest_pixel,loneast_pixel), "Ti")
    indims = InDims("time")
    outdims = OutDims(Dim{:rlevels}(rlevels))
    # outdims = OutDims(outtype=Vector{ReturnLevel{MaximumLikelihoodAbstractExtremeValueModel{ThresholdExceedance}}})

    return mapCube(rlevels_cube, ds, indims=indims, outdims=outdims, threshold=threshold, rlevels=rlevels, minimalvalue=minimalvalue, nthreads=Threads.nthreads())

end

# function moving_fct(cube::YAXArray; fct::Function=mean, latsouth_pixel=1, latnorth_pixel=1, lonwest_pixel=1, loneast_pixel=1)
    
#     indims = InDims(MovingWindow("latitude", latsouth_pixel,latnorth_pixel),
#         MovingWindow("longitude", lonwest_pixel,loneast_pixel))
#     outdims=OutDims()
#     mapCube(moving_fct, cube; fct=fct, indims=indims, outdims=outdims)
# end



function gpfit_cube(xout, xin; threshold=nothing, minimalvalue=1.0)

    # @show size(xin)

    if all(ismissing, xin)
        xout .= missing
        return    
    end

    if isnothing(threshold)
        dataforquantile = xin[xin .> minimalvalue]
        if  !isempty(dataforquantile)
            threshold = quantile(dataforquantile, 0.95)
            exceedances = xin[xin .> threshold]
            xout .= Extremes.gpfit(exceedances)
        else
            threshold = missing
            xout .= missing
            return
        end
        
    end


end


function gpfit_cube(ds::YAXArray; threshold=nothing)

    indims = InDims("time")
    outdims = OutDims(outtype=Union{Missing, MaximumLikelihoodAbstractExtremeValueModel})

    return mapCube(gpfit_cube, ds, indims=indims, outdims=outdims, threshold=threshold, nthreads=Threads.nthreads())

end