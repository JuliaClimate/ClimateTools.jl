using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true"

using Documenter, ClimateTools

makedocs(sitename = "ClimateTools.jl",
    doctest = false,
    format = Documenter.HTML(
    prettyurls = CI,
    ),
    pages = [
    "Home" => "index.md",
    "Installation" => "installation.md",
    "Quick Start" => "quickstart.md",
    "Data and Subsetting" => "datasets.md",
    "Interpolation and Regridding" => "interpolation.md",
    "Bias Correction" => "biascorrection.md",
    "Building Climate Scenarios" => "scenarios.md",
    "Indices and Aggregations" => "indices.md",
    "Examples" => "examples.md",
    "Visualization" => "maps.md",
    "Arithmetic and Ensembles" => "interface.md",
    "Validation and Diagnostics" => "validation.md",
    "Troubleshooting" => "troubleshooting.md",
    "API Overview" => "functions.md",
    "References" => "references.md",
    ]

)


if CI
    deploydocs(
    repo   = "github.com/JuliaClimate/ClimateTools.jl.git",    
    target = "build"
    )
end
