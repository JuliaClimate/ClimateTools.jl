using ClimateTools
using Test
using YAXArrays
using DimensionalData

_xmap_sum_kernel(xout, xin; scale=1.0) = xout .= sum(xin) * scale

function _xmap_pair_kernel(xout, x1, x2, offset; scale=1.0)
    xout .= (sum(x1) + sum(x2) + offset) * scale
end

@testset "xmap utility helpers" begin
    @test ClimateTools._normalize_xmap_dim(:time) == :time
    @test ClimateTools._normalize_xmap_dim("time") == :time
    @test ClimateTools._normalize_xmap_dims(:time) == (:time,)
    @test ClimateTools._normalize_xmap_dims((:longitude, "time")) == (:longitude, :time)
    @test ClimateTools._normalize_xmap_dims(["longitude", :latitude]) == (:longitude, :latitude)

    axis = Dim{:stat}(["sum"])
    @test ClimateTools._normalize_xmap_output_axes() == ()
    @test ClimateTools._normalize_xmap_output_axes(axis) == (axis,)
    @test ClimateTools._normalize_xmap_output_axes((axis,)) == (axis,)
    @test ClimateTools._normalize_xmap_output_axes([axis]) == (axis,)

    cube = YAXArray(
        (Dim{:longitude}([10.0, 20.0]), Dim{:latitude}([45.0, 50.0]), Dim{:time}(1:3)),
        reshape(Float64.(1:12), 2, 2, 3),
    )

    single_window = ClimateTools._xmap_window(cube, (:time,))
    pair_window = ClimateTools._xmap_window(cube, (:longitude, :time))
    @test parentmodule(typeof(single_window)) == YAXArrays.Xmap
    @test parentmodule(typeof(pair_window)) == YAXArrays.Xmap

    spec_no_axes = ClimateTools._xmap_output_spec(())
    spec_with_axes = ClimateTools._xmap_output_spec((axis,); outtype=Float32)
    @test spec_no_axes isa YAXArrays.Xmap.XOutput
    @test spec_with_axes isa YAXArrays.Xmap.XOutput

    dropped = ClimateTools._drop_reduced_xmap_dims(
        YAXArray((Dim{:longitude}([10.0]), Dim{:latitude}([45.0, 50.0])), reshape([1.0, 2.0], 1, 2)),
        (:longitude,),
    )
    @test size(dropped) == (2,)
    @test Tuple(Symbol.(name.(dropped.axes))) == (:latitude,)
end

@testset "xmap call reductions" begin
    cube = YAXArray(
        (Dim{:longitude}([10.0, 20.0]), Dim{:latitude}([45.0, 50.0]), Dim{:time}(1:4)),
        reshape(Float64.(1:16), 2, 2, 4),
    )

    reduced = ClimateTools._xmap_call(
        _xmap_sum_kernel,
        cube;
        reduced_dims="time",
        function_kwargs=(scale=2.0,),
    )
    @test size(reduced) == (2, 2)
    @test Array(reduced) == 2.0 .* dropdims(sum(Array(cube); dims=3); dims=3)

    stat_axis = Dim{:stat}(["combined"])
    shifted = YAXArray(cube.axes, Array(cube) .+ 1.0)
    combined = ClimateTools._xmap_call(
        _xmap_pair_kernel,
        (cube, shifted);
        reduced_dims=[:longitude, "time"],
        output_axes=[stat_axis],
        function_args=(5.0,),
        function_kwargs=(scale=0.5,),
        outtype=Float32,
    )

    # Note: Xmap adds new dimensions (stat) BEFORE the existing ones (latitude)
    expected = Array{Float32}(undef, 1, 2)
    raw_cube = Array(cube)
    raw_shifted = Array(shifted)
    for ilat in 1:2
        expected[1, ilat] = Float32((sum(raw_cube[:, ilat, :]) + sum(raw_shifted[:, ilat, :]) + 5.0) * 0.5)
    end

    @test size(combined) == (1, 2)
    @test Tuple(Symbol.(name.(combined.axes))) == (:stat, :latitude)
    @test Array(combined) ≈ expected
    @test eltype(Array(combined)) == Float32
end
