using ClimateTools
using AxisArrays
# using Lint
using Base.Test

include("indices_test.jl")
include("interface_test.jl")
include("mapping_test.jl")
include("biascorrect_test.jl")
include("indicators_test.jl")

# @test isempty(lintpkg("ClimateTools"))
