function get_dimname(ds::NCDatasets.Dataset, dim::String)

    dimensions = keys(ds.dim)

    found_dim = "NA"

    # try finding with the "axis" attribute

    for idim in dimensions
        if idim != "bnds" && idim != "vertices"
            if haskey(ds[idim].attrib, "axis")
                if ds[idim].attrib["axis"] == dim
                    found_dim = idim
                end
            end
        end
    end

    # if found_dim is still not found, try to find it with the standard_name

    if found_dim == "NA"

        dim_dict = Dict(["T" => "time",
                         "X" => "longitude",
                         "Y" => "latitude"])

        for idim in dimensions
            if idim != "bnds" && idim != "vertices"
                if ds[idim].attrib["standard_name"] == dim_dict[dim]
                    found_dim = idim
                end
            end
        end
    end

    return found_dim
end
