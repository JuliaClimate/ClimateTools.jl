using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"

# PyCall is optional for docs and may be absent from the docs environment.
try
    ENV["PYTHON"] = ""
    Pkg.build("PyCall")
catch
    @info "Skipping PyCall build in docs environment"
end

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
