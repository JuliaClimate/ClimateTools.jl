using Documenter, ClimateTools

makedocs(sitename = "ClimateTools.jl",
    # format = :html,

)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/Balinus/ClimateTools.jl.git",
    julia  = "0.6",
    osname = "linux",
    # deps = nothing,
)
