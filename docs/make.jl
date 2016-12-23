using Documenter, ClimateTools

makedocs(
    doctest = false
)

deploydocs(
    repo = "github.com/Balinus/ClimateTools.jl.git",
    julia = "0.5"
)
