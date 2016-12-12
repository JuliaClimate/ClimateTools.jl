using Documenter, ClimateTools
 
makedocs(modules=[ClimateTools],
        doctest=true)
 
deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/Balinus/ClimateTools.jl.git",
    julia  = "0.5.0",
    osname = "linux")
