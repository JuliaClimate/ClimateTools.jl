"""
    get_dimname(::NCDatasets.Dataset, dim::String)

Scans dataset ds and extract the dimensions X, Y, Z and T. The dataset needs to be CF-compliant. See http://cfconventions.org/.
"""
function get_dimname(ds::NCDatasets.Dataset, dim::String)

    dimensions = keys(ds.dim)

    found_dim = "NA"

    # try finding with the "axis" attribute
    for idim in dimensions
        if idim != "bnds" && idim != "vertices" && idim != "ts" && idim != "string1"
            if haskey(ds[idim].attrib, "axis")
                if ds[idim].attrib["axis"] == dim
                    found_dim = idim
                end
            end
        end
    end

    # if found_dim is still not found, try to find it with the standard_name
    if found_dim == "NA"

        dim_dict = Dict(["T" => ["time"],
                         "Z" => ["pressure", "depth"],
                         "X" => ["longitude"],
                         "Y" => ["latitude"]])

        for idim in dimensions
            if idim != "bnds" && idim != "vertices" && idim != "ts" && idim != "string1"
                attribs = ds[idim].attrib
                if haskey(attribs, "standard_name")
                    name = "standard_name"
                elseif haskey(attribs, "long_name")
                    name = "long_name"
                end
                if in(ds[idim].attrib[name], dim_dict[dim])
                    found_dim = idim
                end
            end
        end
    end

    if found_dim == "NA"
        error(string("Unable to find the dimension ", dim, "!"))
    end

    return found_dim
end

"""
    get_calendar(ds::NCDataset)

Returns the calendar type. See See http://cfconventions.org/.
"""
function get_calendar(ds::NCDataset, timename)
    caltype = ""
    try
        caltype = ds[timename].attrib["calendar"]
    catch
        caltype = "standard"
    end

    return caltype

end
