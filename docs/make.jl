using Documenter, ClimateTools

makedocs(sitename = "ClimateTools.jl",
    pages = [
        "index.md",
        "gettingstarted.md"
        "indices.md",
        "examples.md",
    ]
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/Balinus/ClimateTools.jl.git",
    julia  = "release"
)
