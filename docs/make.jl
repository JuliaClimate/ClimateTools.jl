using Documenter, ClimateTools

makedocs(sitename = "ClimateTools.jl",
    format = :html,
    pages = [
       "index.md",
       "gettingstarted.md",
       "indices.md",
       "interpolation.md",
       "biascorrection.md",
       "maps.md",
       "interface.md",
       "examples.md",
       "functions.md",
       ]

)

deploydocs(
    repo   = "github.com/Balinus/ClimateTools.jl.git",
    # julia  = "1.0",
    # osname = "linux",
    # target = "build",
    # deps = nothing,
    # make = nothing
)
