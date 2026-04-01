@testset "Bias Correct Extremes YAX" begin
    using DataFrames

    dates = collect(DateTime(1961, 1, 1):Day(1):DateTime(1990, 12, 31))
    n = length(dates)

    t = collect(1:n)
    obs_data = Array{Float64}(undef, 2, 2, n)
    for i in 1:2, j in 1:2
        obs_data[i, j, :] .= 2.0 .+ 0.03 .* t .+ 0.2 * i .+ 0.1 * j .+ 4.0 .* abs.(sin.(t ./ (12.0 + i + j)))
    end
    ref_data = obs_data .* 1.10
    fut_data = obs_data .* 1.05
    fut_data[:, :, 100:50:end] .*= 3

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
    tlong = collect(1:nlong)
    fut_long_data = Array{Float64}(undef, 2, 2, nlong)
    for i in 1:2, j in 1:2
        fut_long_data[i, j, :] .= 1.5 .+ 0.02 .* tlong .+ 0.15 * i .+ 0.08 * j .+ 3.5 .* abs.(cos.(tlong ./ (15.0 + i + j)))
    end
    fut_long_data[:, :, 180:120:end] .*= 2.5
    fut_long = YAXArray((Dim{:lon}(1:2), Dim{:lat}(1:2), Dim{:time}(long_dates)), fut_long_data)
    dext_long = biascorrect_extremes(obs, ref, fut_long; detrend=false, gevparams=gevparams)

    expected_len = count(dt -> Dates.monthday(dt) != (2, 29), long_dates)
    @test length(lookup(dext_long, :time)) == expected_len
    @test all(x -> !ismissing(x) && isfinite(x), Array(dext_long))
end