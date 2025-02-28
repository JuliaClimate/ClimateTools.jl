function quantiles(xout,xin; q)
    
    xout .= NaN # Initialisation
    
    if !all(ismissing, xin)        
        xout .= quantile(skipmissing(xin), q)        
    end
end

function quantiles(cube::YAXArray, q=[0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99], dims::String="number")

    # Dimensions
    indims = InDims(dims)        
    outdims = OutDims(Dim{:Quantiles}(q))    

    mapCube(quantiles, cube, indims=indims, outdims=outdims, q=q)

end