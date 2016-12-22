using Documenter, ClimateTools

makedocs(
    sitename = "ClimateTools.jl",
    modules = [ClimateTools]    
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/Balinus/ClimateTools.jl.git",
    julia  = "0.5",
    osname = "linux"
)
