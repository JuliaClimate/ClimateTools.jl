@testset "Bias Correction YAX Edge Cases" begin
    function make_cube(dimsym::Symbol, data, dates)
        return YAXArray((Dim{:longitude}(1:size(data, 1)), Dim{:latitude}(1:size(data, 2)), Dim{dimsym}(dates)), data)
    end

    dates = collect(DateTime(2000, 1, 1):Day(1):DateTime(2001, 12, 31))
    n = length(dates)
    base = reshape(collect(1.0:(2 * 2 * n)), 2, 2, n)

    obs_ti = make_cube(:Ti, base, dates)
    ref_ti = make_cube(:Ti, base .+ 2.0, dates)
    fut_ti = make_cube(:Ti, base .+ 2.0, dates)
    qq_ti = qqmap(obs_ti, ref_ti, fut_ti; method="additive", detrend=false)
    @test :time in name.(qq_ti.axes)
    @test !(:Ti in name.(qq_ti.axes))
    @test length(lookup(qq_ti, :time)) == 730
    @test all(Dates.monthday.(lookup(qq_ti, :time)) .!= Ref((2, 29)))

    obs = make_cube(:time, base, dates)
    ref = make_cube(:time, base .+ 2.0, dates)
    fut = make_cube(:time, base .+ 2.0, dates)

    qq_bulk = qqmap_bulk(obs, ref, fut; method="additive", detrend=false, order=2)
    @test :time in name.(qq_bulk.axes)
    @test length(lookup(qq_bulk, :time)) == 730
    @test all(x -> !ismissing(x) && isfinite(x), Array(qq_bulk))

    short_dates = collect(DateTime(2001, 1, 1):Day(1):DateTime(2001, 12, 31))
    nan_obs_data = reshape(collect(1.0:(1 * 1 * length(short_dates))), 1, 1, length(short_dates))
    nan_ref_data = copy(nan_obs_data)
    nan_fut_data = copy(nan_obs_data)
    nan_obs_data[1, 1, 1:2:end] .= NaN
    nan_ref_data[1, 1, 1:2:end] .= NaN
    nan_fut_data[1, 1, 1:2:end] .= NaN

    nan_obs = make_cube(:time, nan_obs_data, short_dates)
    nan_ref = make_cube(:time, nan_ref_data, short_dates)
    nan_fut = make_cube(:time, nan_fut_data, short_dates)

    qq_nan = qqmap(nan_obs, nan_ref, nan_fut; method="additive", detrend=false, keep_original=false)
    @test all(x -> x isa Number && isnan(x), Array(qq_nan))

    qq_keep = qqmap(nan_obs, nan_ref, nan_fut; method="additive", detrend=false, keep_original=true)
    keep_arr = Array(qq_keep)
    fut_arr = Array(nan_fut)
    @test all(.!ismissing.(keep_arr))
    keep_vec = vec(keep_arr)
    fut_vec = vec(fut_arr)
    @test length(keep_vec) == length(fut_vec)
    @test all(eachindex(keep_vec)) do i
        x = keep_vec[i]
        y = fut_vec[i]
        isnan(y) ? isnan(x) : x == y
    end

    no_time = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2), Dim{:level}(1:3)), rand(2, 2, 3))
    @test_throws ErrorException qqmap(no_time, no_time, no_time)
end