using ClimateTools
using AxisArrays
# using Lint
using Test

include("interface_test.jl")
include("indices_test.jl")
include("mapping_test.jl")
include("biascorrect_test.jl")
include("indicators_test.jl")
include("extract_test.jl")

# @test isempty(lintpkg("ClimateTools"))
