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
    # deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/Balinus/ClimateTools.jl.git",
    julia  = "0.6",
    osname = "linux",
    deps = nothing,
)
