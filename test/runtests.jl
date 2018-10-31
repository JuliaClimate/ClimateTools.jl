using ClimateTools
using AxisArrays
# const axes = Base.axes
using PyPlot
using Shapefile
using NetCDF
using NCDatasets
using Dates
using Statistics
using Random
using ArgCheck
using Unitful: K, °C, m, mm, s, kg
using Unitful: @u_str, ustrip, uconvert
import Unitful
using Test

# println("Running data extraction tests...")
include("extract_test.jl")
# println("Running bias correction tests...")
include("biascorrect_test.jl")
# println("Running interface tests...")
include("interface_test.jl")
# println("Running indices tests...")
include("indices_test.jl")
# include("indicators_test.jl")
# println("Running mapping tests...")
include("mapping_test.jl")
# println("Running functions tests...")
include("functions_test.jl")

# @test isempty(lintpkg("ClimateTools"))
