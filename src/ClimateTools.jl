__precompile__()

module ClimateTools

# External modules
using NetCDF

# Included files
include("functions.jl")
include("indices.jl")

# Exported functions
export windnr, leftorright, inpoly, meshgrid, boxcar3, prcp1, frostdays, summerdays, icingdays, tropicalnights, customthresover, customthresunder, annualmax, annualmin

end #module
