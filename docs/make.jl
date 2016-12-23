#push!(LOAD_PATH,"../src/")

using Documenter, ClimateTools

#makedocs(modules=ClimateTools)

makedocs(
    sitename = "ClimateTools.jl",
    modules = [ClimateTools],
    format = :html,
    clean = false,
    pages = Any["Home" => "index.md"],
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    target = "build",
    repo = "github.com/Balinus/ClimateTools.jl.git",
    julia  = "0.5",
    osname = "linux"
)

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
