#__precompile__()

module ClimateTools

# External modules
using NetCDF
using GMT

# Exported functions
export windnr, leftorright, inpoly, meshgrid, boxcar3, prcp1, frostdays, summerdays, icingdays, tropicalnights, customthresover, customthresunder, annualmax, annualmin

# Included files
include("functions.jl")
include("indices.jl")



end #module
