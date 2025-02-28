function getsymbols(C::ClimGrid)
    latsymbol = Symbol(C.dimension_dict["lat"])
    lonsymbol = Symbol(C.dimension_dict["lon"])

    return latsymbol, lonsymbol
end




