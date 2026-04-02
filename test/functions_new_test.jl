using Test
using ClimateTools
using YAXArrays
using DimensionalData
using Dates
using Statistics

# ---------------------------------------------------------------------------
# Helper: build a simple (lon × lat × time) cube with daily data
# ---------------------------------------------------------------------------
function _make_daily_cube(; nlon=3, nlat=2, start=Date(2020,1,1), stop=Date(2020,3,31))
    lons = collect(1.0:nlon)
    lats = collect(1.0:nlat)
    times = collect(start:Day(1):stop)
    nt = length(times)
    data = Float64[lon + lat + 0.1 * t for lon in lons, lat in lats, t in 1:nt]
    return YAXArray(
        (Dim{:longitude}(lons), Dim{:latitude}(lats), Dim{:time}(times)),
        data,
    )
end

# ---------------------------------------------------------------------------
# Helper: build a sub-daily (hourly) cube for daily_fct / ERA5Land_dailysum
# ---------------------------------------------------------------------------
function _make_hourly_cube(; nlon=2, nlat=2, start=DateTime(2020,1,1,0), ndays=3)
    lons = collect(1.0:nlon)
    lats = collect(1.0:nlat)
    # 24 hourly steps per day
    times = collect(start:Hour(1):start + Day(ndays) - Hour(1))
    nt = length(times)
    data = Float64[lon + lat + sin(2π * t / 24) for lon in lons, lat in lats, t in 1:nt]
    return YAXArray(
        (Dim{:longitude}(lons), Dim{:latitude}(lats), Dim{:time}(times)),
        data,
    )
end


# ===== Date builders =====
@testset "dates_builder_yearmonth" begin
    tuples = [(2020, 1), (2020, 6), (2021, 12)]
    result = dates_builder_yearmonth(tuples)
    @test result == [DateTime(2020, 1), DateTime(2020, 6), DateTime(2021, 12)]
    @test eltype(result) == DateTime
end

@testset "dates_builder_yearmonth_hardcode" begin
    years = [(2020,), (2021,), (2022,)]
    result = dates_builder_yearmonth_hardcode(years, 7)
    @test result == [DateTime(2020, 7), DateTime(2021, 7), DateTime(2022, 7)]
end

@testset "dates_builder_yearmonthday" begin
    tuples = [(2020, 3, 15), (2020, 12, 31)]
    result = dates_builder_yearmonthday(tuples)
    @test result == [DateTime(2020, 3, 15), DateTime(2020, 12, 31)]
end

@testset "dates_builder_yearmonthday_hardcode" begin
    years = [(2020,), (2021,)]
    result = dates_builder_yearmonthday_hardcode(years; imois=6, iday=15)
    @test result == [DateTime(2020, 6), DateTime(2021, 6)]
end


# ===== m2mm =====
@testset "m2mm" begin
    cube = _make_daily_cube(nlon=2, nlat=2, start=Date(2020,1,1), stop=Date(2020,1,5))
    result = m2mm(cube)
    @test Array(result) ≈ Array(cube) .* 1000.0
end


# ===== subsample =====
@testset "subsample" begin
    cube = _make_daily_cube(start=Date(2020,1,1), stop=Date(2020,12,31))

    # Simple range: March to May
    sub = subsample(cube; month1=3, month2=5)
    months_in_sub = month.(collect(sub.time))
    @test all(m -> 3 <= m <= 5, months_in_sub)

    # Year-overlap: November to February
    sub2 = subsample(cube; month1=11, month2=2)
    months_in_sub2 = month.(collect(sub2.time))
    @test all(m -> m >= 11 || m <= 2, months_in_sub2)
    @test length(months_in_sub2) > 0
end


# ===== diff =====
@testset "diff" begin
    times = collect(Date(2020,1,1):Day(1):Date(2020,1,10))
    data = Float64.(1:10)
    cube = YAXArray((Dim{:time}(times),), data)
    result = ClimateTools.diff(cube)
    arr = Array(result)
    # First element is 0.0, rest should be 1.0
    @test arr[1] ≈ 0.0
    @test all(arr[2:end] .≈ 1.0)
end


# ===== cumsum =====
@testset "cumsum" begin
    times = collect(Date(2020,1,1):Day(1):Date(2020,1,5))
    data = Float64[2.0, 3.0, 1.0, 4.0, 5.0]
    cube = YAXArray((Dim{:time}(times),), data)
    result = ClimateTools.cumsum(cube)
    arr = Array(result)
    @test arr ≈ Base.cumsum(data)
end


# ===== daily_fct =====
@testset "daily_fct" begin
    hourly = _make_hourly_cube(nlon=1, nlat=1, ndays=3)
    result = daily_fct(hourly; fct=mean)
    @test length(result.time) == 3

    result_sum = daily_fct(hourly; fct=sum)
    @test length(result_sum.time) == 3

    # Sum should be larger than mean for positive data
    @test Array(result_sum)[1, 1, 1] > Array(result)[1, 1, 1]
end


# ===== yearly_resample =====
@testset "yearly_resample" begin
    cube = _make_daily_cube(start=Date(2019,1,1), stop=Date(2021,12,31))
    result = yearly_resample(cube; fct=mean)
    # Should have 3 years
    @test length(result.time) == 3

    result_max = yearly_resample(cube; fct=maximum)
    result_min = yearly_resample(cube; fct=minimum)
    @test all(Array(result_max) .>= Array(result))
    @test all(Array(result_min) .<= Array(result))
end


# ===== monthly_resample =====
@testset "monthly_resample" begin
    cube = _make_daily_cube(start=Date(2020,1,1), stop=Date(2020,6,30))
    result = monthly_resample(cube; fct=mean)
    # Jan through Jun = 6 months
    @test length(result.time) == 6

    result_sum = monthly_resample(cube; fct=sum)
    # Sum of daily values > mean for months with more than 1 day
    @test Array(result_sum)[1, 1, 1] > Array(result)[1, 1, 1]
end


# ===== ERA5Land_dailysum =====
@testset "ERA5Land_dailysum" begin
    # ERA5 Land accumulations reset daily; diff + sum gives daily totals
    lons = [1.0]
    lats = [1.0]
    # 3 days × 24 hours, cumulative within each day
    times = collect(DateTime(2020,1,1,1):Hour(1):DateTime(2020,1,3,24))
    nt = length(times)
    # Simulate accumulations: values increase within each day
    data = Float64[mod(t - 1, 24) + 1.0 for _ in lons, _ in lats, t in 1:nt]
    cube = YAXArray(
        (Dim{:longitude}(lons), Dim{:latitude}(lats), Dim{:time}(times)),
        data,
    )
    result = ERA5Land_dailysum(cube; keep_1stday=false)
    @test length(result.time) >= 2
    @test all(isfinite, Array(result))

    result_keep = ERA5Land_dailysum(cube; keep_1stday=true)
    @test length(result_keep.time) >= length(result.time)
end


# ===== quantiles =====
@testset "quantiles" begin
    # Create a cube with a "number" dimension (ensemble members)
    members = collect(1:20)
    data = Float64.(members)  # simple increasing data
    cube = YAXArray((Dim{:number}(members),), data)

    q = [0.25, 0.5, 0.75]
    result = quantiles(cube, q, "number")
    arr = Array(result)
    # Median of 1:20 should be ~10.5
    @test abs(arr[2] - 10.5) < 0.5  # 0.5 quantile
    # Q25 < Q50 < Q75
    @test arr[1] < arr[2] < arr[3]
end

@testset "quantiles with spatial dims" begin
    members = collect(1:10)
    lons = [1.0, 2.0]
    data = Float64[m + lon for lon in lons, m in members]
    cube = YAXArray((Dim{:longitude}(lons), Dim{:number}(members)), data)

    q = [0.1, 0.5, 0.9]
    result = quantiles(cube, q, "number")
    @test size(result) == (length(q), length(lons))
    # Higher lon should have higher quantiles
    arr = Array(result)
    @test arr[2, 2] > arr[2, 1]
end


# ===== ensemble_stats =====
@testset "ensemble_stats" begin
    times = collect(Date(2020,1,1):Day(1):Date(2020,1,10))
    data = Float64.(1:10)
    cube = YAXArray((Dim{:time}(times),), data)

    result = ensemble_stats(cube; dim="time")
    arr = Array(result)
    @test arr[1] ≈ 10.0  # max
    @test arr[2] ≈ 5.5   # mean
    @test arr[3] ≈ 1.0   # min
end

@testset "ensemble_stats with spatial dims" begin
    lons = [1.0, 2.0]
    times = collect(Date(2020,1,1):Day(1):Date(2020,1,5))
    data = Float64[lon * t for lon in lons, t in 1:5]
    cube = YAXArray((Dim{:longitude}(lons), Dim{:time}(times)), data)

    result = ensemble_stats(cube; dim="time")
    arr = Array(result)
    # output shape is (stats=3, longitude=2)
    @test size(arr) == (3, 2)
    # lon=1: max=5, mean=3, min=1
    @test arr[1, 1] ≈ 5.0
    @test arr[2, 1] ≈ 3.0
    @test arr[3, 1] ≈ 1.0
    # lon=2: max=10, mean=6, min=2
    @test arr[1, 2] ≈ 10.0
    @test arr[2, 2] ≈ 6.0
    @test arr[3, 2] ≈ 2.0
end

@testset "ensemble_fct alias" begin
    times = collect(Date(2020,1,1):Day(1):Date(2020,1,5))
    data = Float64.(1:5)
    cube = YAXArray((Dim{:time}(times),), data)
    r1 = ensemble_stats(cube)
    r2 = ensemble_fct(cube)
    @test Array(r1) ≈ Array(r2)
end


# ===== autocorrelation =====
@testset "autocorrelation" begin
    # Create a simple time series
    n = 200
    times = collect(Date(2000,1,1):Day(1):Date(2000,1,1) + Day(n-1))
    data = Float64[sin(2π * t / 30) + 0.1 * t for t in 1:n]
    cube = YAXArray((Dim{:time}(times),), data)

    result = autocorrelation(cube; lags=10)
    arr = Array(result)
    @test length(arr) == 10
    # Lag-1 autocorrelation of a smooth signal should be positive and high
    @test arr[1] > 0.5
end


# ===== hurst =====
@testset "hurst" begin
    n = 500
    times = collect(Date(2000,1,1):Day(1):Date(2000,1,1) + Day(n-1))
    data = [sin(2 * π * t / 30) + 0.1 * cos(2 * π * t / 7) for t in 1:n]
    cube = YAXArray((Dim{:time}(times),), data)

    result = hurst(cube; k=10)
    h = Array(result)[]
    @test isfinite(h)
end


# ===== rlevels_cube =====
@testset "rlevels_cube" begin
    n = 3650  # ~10 years of daily data
    times = collect(Date(2000,1,1):Day(1):Date(2000,1,1) + Day(n-1))
    # Generate positive precipitation-like data (exponential-ish)
    rng = 42  # seed-like offset
    data = Float64[max(0.0, 5.0 * sin(2π * t / 365) + 3.0 + 0.5 * ((t * 7 + rng) % 17 - 8)) for t in 1:n]
    cube = YAXArray((Dim{:time}(times),), data)

    rlevels = [2, 10, 50]
    result = rlevels_cube(cube; rlevels=rlevels, minimalvalue=1.0)
    arr = Array(result)
    @test length(arr) == 3
    # Return levels should increase with return period
    @test all(isfinite, arr)
    @test arr[1] <= arr[2] <= arr[3]
end

@testset "rlevels_cube all missing" begin
    n = 100
    times = collect(Date(2000,1,1):Day(1):Date(2000,1,1) + Day(n-1))
    data = fill(missing, n)
    cube = YAXArray((Dim{:time}(times),), data)

    rlevels = [2, 10]
    result = rlevels_cube(cube; rlevels=rlevels)
    arr = Array(result)
    @test all(ismissing, arr)
end
