# """
#     getdim_lat(ds::NCDatasets.Dataset)

# Returns the name of the "latitude" dimension and the status related to a regular grid. The latitude dimension is usually "latitude", "lat", "y", "yc", "rlat".
# """
# function getdim_lat(ds::NCDatasets.Dataset)

#     if sum(keys(ds.dim) .== "rlat") == 1
#         return "rlat", true
#     elseif sum(keys(ds.dim) .== "lat") == 1
#         return "lat", false
#     elseif sum(keys(ds.dim) .== "latitude") == 1
#         return "latitude", false
#     elseif sum(keys(ds.dim) .== "y") == 1
#         return "y", true
#     elseif sum(keys(ds.dim) .== "yc") == 1
#         return "yc", true
#     elseif sum(keys(ds.dim) .== "lat_c") == 1
#         return "lat_c", false
#     else
#         error("Manually verify x/lat dimension name")
#     end

# end

# """
#     getdim_lon(ds::NCDatasets.Dataset)

# Returns the name of the "longitude" dimension and the status related to a regular grid. The longitude dimension is usually "longitue", "lon", "x", "xc", "rlon".
# """
# function getdim_lon(ds::NCDatasets.Dataset)

#     if sum(keys(ds.dim) .== "rlon") == 1
#         return "rlon", true
#     elseif sum(keys(ds.dim) .== "lon") == 1
#         return "lon", false
#     elseif sum(keys(ds.dim) .== "longitude") == 1
#         return "longitude", false
#     elseif sum(keys(ds.dim) .== "x") == 1
#         return "x", false
#     elseif sum(keys(ds.dim) .== "xc") == 1
#         return "xc", false
#     elseif sum(keys(ds.dim) .== "lon_c") == 1
#         return "lon_c", false
#     else
#         error("Manually verify x/lat dimension name")
#     end

# end

# """
#     latgridname(ds::NCDatasets.Dataset)

# Returns the name of the latitude grid when datasets is not on a rectangular grid.
# """
# function latgridname(ds::NCDatasets.Dataset)

#     # if in("lat", keys(ds))
#     #     return "lat"
#     # elseif in("latitude", keys(ds))
#     #     return "latitude"
#     # elseif in("lat_c", keys(ds))
#     #     return "lat_c"
#     # else
#     #     error("Variable name is not supported. File an issue on https://github.com/Balinus/ClimateTools.jl/issues")
#     # end

#     names = ["degree_north",
#              "degrees_north",
#              "degreeN",
#              "degreesN",
#              "degree_N",
#              "degrees_N"]

#     varnames = keys(ds)

#     found_var = "NA"

#     for ivar in varnames

#         if isdefined(ds[ivar], :attrib)
#             if haskey(ds[ivar].attrib, "units")
#                 if in(ds[ivar].attrib["units"], names)
#                     found_var = ivar
#                 end
#             end
#         end

#     end

#     return found_var

# end

# """
#     longridname(ds::NCDatasets.Dataset)

# Returns the name of the longitude grid when datasets is not on a rectangular grid.
# """
# function longridname(ds::NCDatasets.Dataset)

#     # if in("lon", keys(ds))
#     #     return "lon"
#     # elseif in("longitude", keys(ds))
#     #     return "longitude"
#     # elseif in("lon_c", keys(ds))
#     #     return "lon_c"
#     # else
#     #     error("Variable name is not supported. File an issue on https://github.com/Balinus/ClimateTools.jl/issues")
#     # end

#     names = ["degree_east",
#              "degrees_east",
#              "degreeE",
#              "degreesE",
#              "degree_E",
#              "degrees_E"]

#     varnames = keys(ds)

#     found_var = "NA"

#     for ivar in varnames

#         if isdefined(ds[ivar], :attrib)
#             if haskey(ds[ivar].attrib, "units")
#                 if in(ds[ivar].attrib["units"], names)
#                     found_var = ivar
#                 end
#             end
#         end

#     end

#     return found_var



# end

# """
#     getdim_tim(ds::NCDatasets.Dataset)

# Returns the name of the "time" dimension.
# """
# function getdim_tim(ds::NCDatasets.Dataset)

#     if sum(keys(ds.dim) .== "time") == 1
#         return "time", false
#     elseif sum(keys(ds.dim) .== "tim") == 1
#         return "tim", false
#     else
#         error("Manually verify time dimension name")
#     end

# end

# """
#     extractdata(data, msk, idxtimebeg, idxtimeend)

# Returns the data contained in netCDF file, using the appropriate mask and time index. Used internally by `load`.

# """
# function extractdata(data, msk, idxtimebeg, idxtimeend)

#     # idlon, idlat = findn(.!isnan.(msk))
#     begin
#         I = Base.findall(!isnan, msk)
#         idlon, idlat = (getindex.(I, 1), getindex.(I, 2))
#     end
#     minXgrid = minimum(idlon)
#     maxXgrid = maximum(idlon)
#     minYgrid = minimum(idlat)
#     maxYgrid = maximum(idlat)

#     if ndims(data) == 3
#         dataout = data[minXgrid:maxXgrid, minYgrid:maxYgrid, idxtimebeg:idxtimeend]
#         # Permute dims
#         # data = permutedims(data, [3, 1, 2])
#     elseif ndims(data) == 4
#         dataout = data[minXgrid:maxXgrid, minYgrid:maxYgrid, :, idxtimebeg:idxtimeend]
#         # Permute dims
#         # data = permutedims(data, [4, 1, 2, 3])
#     end

#     return dataout
# end

# function extractdata2D(data, msk)

#     # idlon, idlat = findn(.!isnan.(msk))
#     begin
#         I = Base.findall(!isnan, msk)
#         idlon, idlat = (getindex.(I, 1), getindex.(I, 2))
#     end
#     minXgrid = minimum(idlon)
#     maxXgrid = maximum(idlon)
#     minYgrid = minimum(idlat)
#     maxYgrid = maximum(idlat)

#     data = data[minXgrid:maxXgrid, minYgrid:maxYgrid]

#     return data

# end


# """
#     get_mapping(ds::Array{String,1})

# Returns the grid_mapping of Dataset *ds*
# """
# function get_mapping(K::Array{String,1})

#     if in("rotated_pole", K)
#         return "rotated_pole"

#     elseif in("lambert_conformal_conic", K)
#         return "lambert_conformal_conic"

#     elseif in("rotated_latitude_longitude", K)
#         return "rotated_latitude_longitude"

#     elseif in("rotated_mercator", K)
#         return "rotated_mercator"

#     elseif in("crs", K)
#         return "crs"

#     elseif in("polar_stereographic", K)
#         return "polar_stereographic"

#     else
#         return "Regular_longitude_latitude"
#     end

# end

# """
#     get_mapping(ds::Array{String,1})

# Returns the grid_mapping of Dataset *ds*
# """
# function get_mapping(ds::NCDatasets.Dataset, vari)

#     grid = ""

#     try
#         grid = ds[vari].attrib["grid_mapping"]
#     catch
#         grid = "Regular_longitude_latitude"
#     end

#     return grid


# end

# function build_grid_mapping(ds::NCDatasets.Dataset, grid_mapping::String)

#     if ClimateTools.@isdefined grid_mapping
#         # map_dim = varattrib[grid_mapping]
#         map_attrib = Dict(ds[grid_mapping].attrib)
#         map_attrib["grid_mapping"] = grid_mapping
#     else
#         error("File an issue on https://github.com/Balinus/ClimateTools.jl/issues to get the grid supported")
#     end

#     return map_attrib

# end
