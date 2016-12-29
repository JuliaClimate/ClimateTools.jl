__precompile__()

module ClimateTools

# External modules
using NetCDF

# Exported functions
export windnr, leftorright, inpoly, meshgrid, boxcar3, prcp1, frostdays, summerdays, icingdays, tropicalnights, customthresover, customthresunder, annualmax, annualmin, netcdf2julia

# Included files
include("functions.jl")
include("indices.jl")
include("extract.jl")


end #module
