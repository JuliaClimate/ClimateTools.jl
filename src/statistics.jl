function quantiles(xout,xin; q)
    
    xout .= NaN # Initialisation
    
    if !all(ismissing, xin)        
        xout .= quantile(skipmissing(xin), q)        
    end
end

function quantiles(cube::YAXArray, q=[0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99], dims::String="number")
    return _xmap_call(
        quantiles,
        cube;
        reduced_dims=dims,
        output_axes=(Dim{:Quantiles}(q),),
        function_kwargs=(q=q,),
    )

end