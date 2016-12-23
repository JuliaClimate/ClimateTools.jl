push!(LOAD_PATH,"../src/")

using Documenter, ClimateTools

makedocs(modules = ClimateTools)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/Balinus/ClimateTools.jl.git",
    julia  = "0.5",
    osname = "linux"
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
