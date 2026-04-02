@testset "Bias Correct Extremes YAX" begin
    using DataFrames

    dates = collect(DateTime(1961, 1, 1):Day(1):DateTime(1990, 12, 31))
    n = length(dates)

    rainfall_template = [
        0.0, 2.3, 1.3, 6.9, 4.6, 0.0, 1.0, 1.5, 1.8, 1.8, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        3.3, 2.8, 3.0, 0.0, 8.1, 1.5, 4.1, 0.0, 0.0, 0.0, 0.0, 0.0,
        4.8, 31.8, 0.0, 1.5, 25.4, 5.1, 15.0, 16.8, 16.3, 0.0, 0.0, 11.7,
        2.3, 2.0, 10.9, 8.1, 2.3, 1.5, 0.0, 0.0, 0.0, 3.0, 1.8, 2.5,
        3.0, 6.6, 2.0, 8.4, 7.4, 11.9, 32.5, 10.7, 2.5, 18.3, 5.1, 13.5,
        10.9, 8.1, 0.8, 12.7, 0.3, 15.7, 18.5, 2.3, 1.8, 5.8, 2.0, 7.1,
        2.3, 0.0, 10.7, 6.9, 4.8, 0.0, 3.8, 3.8, 5.8, 8.4, 7.6, 3.0,
        3.6, 3.6, 4.8, 14.7, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 3.0, 14.2, 5.6, 2.5, 1.8, 6.4, 0.8, 2.0, 5.3, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.3, 0.0,
        0.0, 0.0, 0.0, 2.3, 5.1, 0.0, 1.3, 4.6, 0.0, 0.0, 0.0, 0.0,
        0.3, 0.8, 6.4, 17.0, 0.5, 5.1, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 15.5, 1.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 3.0, 3.8, 0.0, 17.8, 0.0, 13.0, 8.1, 0.0, 0.0, 0.0, 0.0,
        0.0, 5.1, 0.0, 2.3, 1.3, 0.5, 0.3, 27.9, 0.0, 0.0, 0.0, 4.6,
        2.0, 3.0, 7.9, 1.0, 0.0, 4.6, 0.0, 20.3, 14.7, 1.0, 7.6, 3.6,
        0.0, 3.3, 7.1, 4.1, 2.5, 0.0, 0.0, 0.0, 0.0, 24.1, 4.3, 0.0,
        0.0, 0.0, 0.0, 4.1, 0.0, 0.5, 0.0, 3.6, 17.5, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.6, 0.0, 0.0, 7.4, 0.0, 11.7,
        7.6, 1.0, 6.6, 11.4, 1.0, 3.6, 6.6, 1.3, 7.6, 1.8, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.8, 1.3, 0.8,
        0.0, 0.0, 0.0, 0.0, 0.0, 8.1, 15.0, 2.0, 0.3, 11.9, 1.3, 5.1,
        10.2, 1.8, 3.3, 14.7, 13.0, 0.0, 0.0, 9.4, 0.0, 2.5, 0.0, 0.0,
        0.0, 0.0, 2.0, 3.3, 11.4, 13.2, 7.9, 0.0, 0.0, 0.0, 4.6, 0.0,
        0.0, 0.0, 0.0, 2.5, 3.0, 8.9, 6.4, 6.4, 22.9, 15.7, 0.0, 16.8,
        3.8, 22.1, 5.1, 11.9, 6.9, 0.0, 0.0, 0.0, 11.9, 2.3, 4.3, 4.6,
        18.3, 0.0, 31.8, 8.1, 1.5, 5.8, 3.0, 0.0, 0.0, 0.0, 11.4, 21.3,
        11.7, 20.3, 1.3, 44.5, 14.0,
    ]
    obs_template = repeat(rainfall_template, cld(n, length(rainfall_template)))[1:n]

    obs_data = Array{Float64}(undef, 2, 2, n)
    for i in 1:2, j in 1:2
        scale = 1.0 + 0.04 * (i - 1) + 0.03 * (j - 1)
        obs_data[i, j, :] .= obs_template .* scale
    end
    ref_data = obs_data .* 1.05
    fut_data = obs_data .* 1.08
    fut_data[:, :, 200:365:n] .*= 1.2

    obs = YAXArray((Dim{:lon}(1:2), Dim{:lat}(1:2), Dim{:time}(dates)), obs_data)
    ref = YAXArray((Dim{:lon}(1:2), Dim{:lat}(1:2), Dim{:time}(dates)), ref_data)
    fut = YAXArray((Dim{:lon}(1:2), Dim{:lat}(1:2), Dim{:time}(dates)), fut_data)

    qq = qqmap(obs, ref, fut; method="multiplicative", detrend=false)
    dext = biascorrect_extremes(obs, ref, fut; detrend=false)

    @test size(dext) == size(qq)
    @test :time in name.(dext.axes)
    @test all(x -> !ismissing(x) && isfinite(x), Array(dext))
    @test count(abs.(Array(dext) .- Array(qq)) .> 1e-8) > 0

    gevparams = DataFrame(
        lat=[1.0, 1.0, 2.0, 2.0],
        lon=[1.0, 2.0, 1.0, 2.0],
        mu=[20.0, 22.0, 21.0, 23.0],
        sigma=[5.0, 4.0, 5.5, 4.5],
        xi=[0.10, 0.05, 0.08, 0.03],
    )

    dext_gev = biascorrect_extremes(obs, ref, fut; detrend=false, gevparams=gevparams)
    @test size(dext_gev) == size(qq)
    @test all(x -> !ismissing(x) && isfinite(x), Array(dext_gev))
    @test count(abs.(Array(dext_gev) .- Array(qq)) .> 1e-8) > 0

    dext_detrend = biascorrect_extremes(obs, ref, fut; detrend=true, gevparams=gevparams)
    @test size(dext_detrend) == size(qq)

    long_dates = collect(DateTime(1961, 1, 1):Day(1):DateTime(2020, 12, 30))
    nlong = length(long_dates)
    fut_long_template = repeat(rainfall_template, cld(nlong, length(rainfall_template)))[1:nlong]
    fut_long_data = Array{Float64}(undef, 2, 2, nlong)
    for i in 1:2, j in 1:2
        scale = (1.0 + 0.04 * (i - 1) + 0.03 * (j - 1)) * 1.06
        fut_long_data[i, j, :] .= fut_long_template .* scale
    end
    fut_long_data[:, :, 220:400:nlong] .*= 1.15
    fut_long = YAXArray((Dim{:lon}(1:2), Dim{:lat}(1:2), Dim{:time}(long_dates)), fut_long_data)
    dext_long = biascorrect_extremes(obs, ref, fut_long; detrend=false, gevparams=gevparams)

    expected_len = count(dt -> Dates.monthday(dt) != (2, 29), long_dates)
    @test length(lookup(dext_long, :time)) == expected_len
    @test all(x -> !ismissing(x) && isfinite(x), Array(dext_long))
end