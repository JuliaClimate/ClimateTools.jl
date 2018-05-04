using Documenter, ClimateTools

makedocs(sitename = "ClimateTools.jl",
    format = :html,
    pages = [
        "Home" => "index.md",
        "Getting Started" => "gettingstarted.md",
        "Climate Indices" => "indices.md",
        "Interpolation" => "interpolation.md",
        "Bias correction" => "biascorrection.md",
        "Maps" => "maps.md",
        "Interface" => "interface.md",
        "Examples" => "examples.md",
    ]
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/Balinus/ClimateTools.jl.git",
    julia  = "0.6",
    osname = "linux",
)
