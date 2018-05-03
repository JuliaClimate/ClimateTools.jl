using Documenter, ClimateTools

makedocs()

deploydocs(
   deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/Balinus/ClimateTools.jl.git",
    target = "build",
    # deps = nothing,
    make = nothing,
    julia  = "0.6",
#    osname = "linux",
)
