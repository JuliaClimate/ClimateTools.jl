__precompile__()

module ClimateTools
using NetCDF
include("functions.jl")
include("indices.jl")

export windnr, leftorright, inpoly, meshgrid, boxcar3, prcp1, frostdays, summerdays, icingdays, tropicalnights

end #module
