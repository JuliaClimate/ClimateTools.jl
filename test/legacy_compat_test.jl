@testset "Legacy Compatibility" begin
    function first_gridpoint_timeseries(c)
        axes_names = name.(c.axes)
        lon_idx = findfirst(==(:longitude), axes_names)
        lat_idx = findfirst(==(:latitude), axes_names)
        arr = Array(c)
        idx = ntuple(i -> (i == lon_idx || i == lat_idx) ? 1 : Colon(), ndims(arr))
        return vec(arr[idx...])
    end

    dates = collect(DateTime(2000, 1, 1):Day(1):DateTime(2002, 12, 31))

    yearly_data = Array{Float64}(undef, 2, 2, length(dates))
    yearly_data[1,1,:] = vcat(fill(1.0, 366), fill(2.0, 365), fill(3.0, 365))
    yearly_data[1,2,:] = yearly_data[1,1,:]
    yearly_data[2,1,:] = yearly_data[1,1,:]
    yearly_data[2,2,:] = yearly_data[1,1,:]

    cube = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2), Dim{:time}(dates)), yearly_data)

    @test first_gridpoint_timeseries(annualmax(cube)) == [1.0, 2.0, 3.0]
    @test first_gridpoint_timeseries(annualmin(cube)) == [1.0, 2.0, 3.0]
    @test first_gridpoint_timeseries(annualmean(cube)) == [1.0, 2.0, 3.0]
    @test first_gridpoint_timeseries(annualsum(cube)) == [366.0, 730.0, 1095.0]

    threshold_cube = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2), Dim{:time}(dates)), yearly_data .* 10)
    @test first_gridpoint_timeseries(summerdays(threshold_cube)) == [0.0, 0.0, 365.0]
    @test first_gridpoint_timeseries(tropicalnights(threshold_cube)) == [0.0, 0.0, 365.0]
    @test first_gridpoint_timeseries(customthresover(threshold_cube, 15)) == [0.0, 365.0, 365.0]
    @test first_gridpoint_timeseries(customthresunder(threshold_cube, 25)) == [366.0, 365.0, 0.0]

    precip_cube = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2), Dim{:time}(dates)), yearly_data)
    @test first_gridpoint_timeseries(prcp1(precip_cube)) == [366.0, 365.0, 365.0]

    cold_cube = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2), Dim{:time}(dates)), -yearly_data)
    @test first_gridpoint_timeseries(frostdays(cold_cube)) == [366.0, 365.0, 365.0]
    @test first_gridpoint_timeseries(icingdays(cold_cube)) == [366.0, 365.0, 365.0]

    six_hour_dates = collect(DateTime(2000, 1, 1):Hour(6):DateTime(2000, 1, 3, 18))
    six_hour_data = reshape(collect(1.0:length(six_hour_dates)), 1, 1, :)
    subdaily = YAXArray((Dim{:longitude}(1:1), Dim{:latitude}(1:1), Dim{:time}(six_hour_dates)), six_hour_data)

    @test first_gridpoint_timeseries(daymean(subdaily)) == [2.5, 6.5, 10.5]
    @test first_gridpoint_timeseries(daysum(subdaily)) == [10.0, 26.0, 42.0]

    indicator_dates = collect(DateTime(2003, 1, 1):Day(1):DateTime(2003, 1, 3))
    tasmin_data = Array{Float64}(undef, 2, 2, 3)
    tasmax_data = Array{Float64}(undef, 2, 2, 3)
    tasmin_data[:,:,1] .= -10.0
    tasmin_data[:,:,2] .= 0.0
    tasmin_data[:,:,3] .= 10.0
    tasmax_data[1,1,:] .= 0.0
    tasmax_data[1,2,:] .= 10.0
    tasmax_data[2,1,:] .= 20.0
    tasmax_data[2,2,:] .= 30.0

    tasmin = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2), Dim{:time}(indicator_dates)), tasmin_data)
    tasmax = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2), Dim{:time}(indicator_dates)), tasmax_data)

    @test Array(meantemperature(tasmin, tasmax))[1,1,:] == [-5.0, 0.0, 5.0]
    @test Array(diurnaltemperature(tasmin, tasmax, 0.4))[2,2,:] == [14.0, 18.0, 22.0]

    huss_data = Array{Float64}(undef, 2, 2, 3)
    huss_data[1,1,:] .= 0.0
    huss_data[1,2,:] .= 0.1
    huss_data[2,1,:] .= 0.5
    huss_data[2,2,:] .= 1.0
    ps_data = Array{Float64}(undef, 2, 2, 3)
    ps_data[:,:,1] .= 10.0
    ps_data[:,:,2] .= 100.0
    ps_data[:,:,3] .= 1000.0

    huss = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2), Dim{:time}(indicator_dates)), huss_data)
    ps = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2), Dim{:time}(indicator_dates)), ps_data)

    vp = vaporpressure(huss, ps)
    @test round(Array(vp)[2,2,3], digits=10) == round(1000.0 / 1.622, digits=10)

    psl = ps
    orog = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2)), [0.0 10.0; 0.0 10.0])
    tas = YAXArray((Dim{:longitude}(1:2), Dim{:latitude}(1:2), Dim{:time}(indicator_dates)), fill(300.0, 2, 2, 3))

    sp = approx_surfacepressure(psl, orog, tas)
    @test size(sp) == size(psl)
    @test Array(sp)[1,1,2] == Array(psl)[1,1,2]

    vp2 = vaporpressure(huss, psl, orog, tas)
    @test size(vp2) == size(psl)

    wbgt_out = wbgt(tas, vp2)
    @test size(wbgt_out) == size(tas)
    @test all(isfinite, Array(wbgt_out))

    @test ClimateTools.timeresolution([1, 2, 3]) == "24h"
    @test ClimateTools.daymean_factor("6h") == 4
    @test ClimateTools.pr_timefactor("3h") == 10800.0
    @test size(ensemble_fct(cube)) == (3, 2, 2)
end