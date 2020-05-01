using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"
ENV["PYTHON"] = ""
Pkg.build("PyCall")
using Documenter, ClimateTools

makedocs(sitename = "ClimateTools.jl",
    doctest = false,
    format = Documenter.HTML(
    prettyurls = CI,
    ),
    pages = [
       "index.md",
       "installation.md",
       "datasets.md",
       "indices.md",
       "interpolation.md",
       "biascorrection.md",
       "maps.md",
       "interface.md",
       "examples.md",
       "functions.md",
       ]

)


if CI
    deploydocs(
    repo   = "github.com/JuliaClimate/ClimateTools.jl.git",    
    target = "build"
    )
end
