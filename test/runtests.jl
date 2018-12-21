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
using ClimateTools


include("interface_test.jl")
include("extract_test.jl")
include("biascorrect_test.jl")
include("indices_test.jl")
include("mapping_test.jl")
include("functions_test.jl")
include("spatial_test.jl")
