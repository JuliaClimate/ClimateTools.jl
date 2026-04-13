@testset "Aggregate" begin
    x = [NaN, 34.5, 23.1, missing, NaN, 56.9, 21.1]
    index_list = [1:4, 5:7]

    out = fill(NaN, length(index_list))
    ClimateTools.daily_fct(out, x; fct=mean, index_list=index_list)
    @test isapprox(out[1], 28.8; atol=1e-10)
    @test isapprox(out[2], 39.0; atol=1e-10)

    out2 = fill(NaN, length(index_list))
    ClimateTools.yearly_resample(out2, x; fct=mean, index_list=index_list)
    @test isapprox(out2[1], 28.8; atol=1e-10)
    @test isapprox(out2[2], 39.0; atol=1e-10)

    hourly_times = collect(DateTime(2020, 1, 1, 0):Hour(1):DateTime(2020, 1, 2, 23))
    hourly_values = Union{Missing, Float64}[
        1, 2, missing, 4, NaN, 6, 7, 8, 9, 10, 11, 12,
        13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
        25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36,
        37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
    ]
    hourly_cube = YAXArray((Dim{:time}(hourly_times),), hourly_values)

    dmmm = ClimateTools.daily_max_mean_min(hourly_cube)
    dmmm_array = Array(dmmm)
    first_day = Float64[1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
    second_day = Float64[25:48;]

    @test size(dmmm) == (2, 3)
    @test Tuple(Symbol.(name.(dmmm.axes))) == (:time, :stats)
    @test collect(lookup(dmmm, :stats)) == ["max", "mean", "min"]
    @test dmmm_array[1, 1] == Base.maximum(first_day)
    @test dmmm_array[1, 2] ≈ mean(first_day)
    @test dmmm_array[1, 3] == Base.minimum(first_day)
    @test dmmm_array[2, 1] == Base.maximum(second_day)
    @test dmmm_array[2, 2] ≈ mean(second_day)
    @test dmmm_array[2, 3] == Base.minimum(second_day)

    quantile_levels = [0.1, 0.5, 0.9]
    time_axis = [Date(2020, 1, 1), Date(2020, 1, 2)]
    lon_axis = [10.0, 20.0]
    lat_axis = [45.0]
    members = collect(1:5)
    percentile_data = Array{Float64}(undef, length(time_axis), length(lon_axis), length(lat_axis), length(members))
    for time_index in eachindex(time_axis), lon_index in eachindex(lon_axis), lat_index in eachindex(lat_axis), member_index in eachindex(members)
        percentile_data[time_index, lon_index, lat_index, member_index] = 100 * time_index + 10 * lon_index + members[member_index]
    end
    percentile_cube = YAXArray(
        (Dim{:time}(time_axis), Dim{:longitude}(lon_axis), Dim{:latitude}(lat_axis), Dim{:number}(members)),
        percentile_data,
    )

    percentile_result = ClimateTools.percentiles(percentile_cube, quantile_levels)
    percentile_array = Array(percentile_result)
    expected_quantiles = Statistics.quantile(collect(111.0:115.0), quantile_levels)

    @test size(percentile_result) == (3, 2, 2, 1)
    @test Tuple(Symbol.(name.(percentile_result.axes))) == (:Quantiles, :time, :longitude, :latitude)
    @test collect(lookup(percentile_result, :Quantiles)) == quantile_levels
    @test percentile_array[:, 1, 1, 1] ≈ expected_quantiles
    @test percentile_array[1, 2, 2, 1] < percentile_array[2, 2, 2, 1] < percentile_array[3, 2, 2, 1]

    six_hour_times = collect(DateTime(2020, 1, 1, 0):Hour(6):DateTime(2020, 1, 3, 18))
    six_hour_values = Float64.(1:length(six_hour_times))
    six_hour_cube = YAXArray((Dim{:time}(six_hour_times),), six_hour_values)

    daily_max_result = ClimateTools.daily_max(six_hour_cube)
    @test size(daily_max_result) == (3,)
    @test Array(daily_max_result) == [4.0, 8.0, 12.0]

    shifted_daily_sum = ClimateTools.daily_fct(six_hour_cube; fct=sum, shifthour=6)
    @test collect(lookup(shifted_daily_sum, :time)) == [
        DateTime(2020, 1, 1),
        DateTime(2020, 1, 2),
        DateTime(2020, 1, 3),
        DateTime(2020, 1, 4),
    ]
    @test Array(shifted_daily_sum) == [6.0, 22.0, 38.0, 12.0]

    yearly_dates = collect(Date(2019, 1, 1):Day(1):Date(2021, 12, 31))
    yearly_cube = YAXArray((Dim{:time}(yearly_dates),), ones(Float64, length(yearly_dates)))

    annual_max_cube = ClimateTools.annualmax(yearly_cube)
    annual_sum_cube = ClimateTools.annualsum(yearly_cube)

    @test collect(lookup(annual_sum_cube, :time)) == [
        DateTime(2019, 7, 1),
        DateTime(2020, 7, 1),
        DateTime(2021, 7, 1),
    ]
    @test Array(annual_max_cube) == [1.0, 1.0, 1.0]
    @test Array(annual_sum_cube) == [365.0, 366.0, 365.0]
end