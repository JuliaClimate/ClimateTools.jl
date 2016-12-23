push!(LOAD_PATH,"../src/")

using Documenter, ClimateTools

makedocs(modules = ClimateTools)

deploydocs(
    repo = "github.com/Balinus/ClimateTools.jl.git"
)

#makedocs(
#    sitename = "ClimateTools.jl",
#    modules = [ClimateTools],
#    format = :html,
#    clean = false,
#    pages = Any["Home" => "index.md"],
#)

#deploydocs(
#    repo   = "github.com/Balinus/ClimateTools.jl.git",
#    target = "build",
#    deps   = nothing,
#    make   = nothing
#)

#deploydocs(
#    target = "build",
#    deps = nothing,
#    make = nothing,
#    repo = "github.com/Balinus/ClimateTools.jl.git",
#)
