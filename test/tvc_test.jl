@testset "Time Variability Correction" begin
    scales = [30, 15, 7, 3]
    n = 240
    time_index = collect(1:n)

    obs = 20 .+ 3 .* sin.(2 .* pi .* time_index ./ 40) .+ 1.5 .* sin.(2 .* pi .* time_index ./ 9)
    raw_train = 18 .+ 1.8 .* sin.(2 .* pi .* time_index ./ 40 .+ 0.2) .+ 0.6 .* sin.(2 .* pi .* time_index ./ 9 .+ 0.1)
    raw_val = 18.5 .+ 1.7 .* sin.(2 .* pi .* (time_index .+ 11) ./ 40 .+ 0.15) .+ 0.55 .* sin.(2 .* pi .* (time_index .+ 11) ./ 9)

    model = fit_tvc(obs, raw_train; scales=scales)
    applied = apply_tvc(model, raw_val)
    direct = tvc(obs, raw_train, raw_val; scales=scales)
    tail_range = (model.init_index + 1):length(applied)

    @test all(isnan, direct[1:model.init_index])
    @test applied[tail_range] ≈ direct[tail_range]
    @test all(isnan, applied[1:model.init_index])
    @test all(isfinite, applied[tail_range])

    obs_matrix = ClimateTools._tvc_truncate(ClimateTools._tvc_decompose(obs, scales), model.init_index)
    raw_matrix = ClimateTools._tvc_truncate(ClimateTools._tvc_decompose(raw_val, scales), model.init_index)
    corrected_matrix = ClimateTools._tvc_apply_components(model, raw_matrix)
    _, obs_covariance = ClimateTools._tvc_stats(obs_matrix)
    _, raw_covariance = ClimateTools._tvc_stats(raw_matrix)
    _, corrected_covariance = ClimateTools._tvc_stats(corrected_matrix)

    @test sum(abs2, corrected_covariance .- obs_covariance) < sum(abs2, raw_covariance .- obs_covariance)

    singular_obs = fill(3.0, 180)
    singular_train = fill(1.0, 180)
    singular_val = fill(1.0, 190)
    singular_model = fit_tvc(singular_obs, singular_train; scales=[20, 10, 5, 2], eig_floor=1e-6)
    singular_corrected = apply_tvc(singular_model, singular_val)
    singular_tail = (singular_model.init_index + 1):length(singular_corrected)

    @test all(isnan, singular_corrected[1:singular_model.init_index])
    @test all(singular_corrected[singular_tail] .≈ 3.0)
end

@testset "Time Variability Correction YAXArray" begin
    function make_cube(dimsym::Symbol, values, dates)
        data = reshape(values, 1, 1, length(values))
        return YAXArray((Dim{:longitude}(1:1), Dim{:latitude}(1:1), Dim{dimsym}(dates)), data)
    end

    scales = [20, 10, 5, 2]
    dates = collect(DateTime(2000, 1, 1):Day(1):DateTime(2001, 12, 31))
    n = length(dates)
    time_index = collect(1:n)

    obs_values = 15 .+ 2.5 .* sin.(2 .* pi .* time_index ./ 30) .+ 1.1 .* cos.(2 .* pi .* time_index ./ 11)
    train_values = 13 .+ 1.9 .* sin.(2 .* pi .* time_index ./ 30 .+ 0.15) .+ 0.7 .* cos.(2 .* pi .* time_index ./ 11 .+ 0.1)
    val_values = 13.5 .+ 1.7 .* sin.(2 .* pi .* (time_index .+ 17) ./ 30 .+ 0.15) .+ 0.8 .* cos.(2 .* pi .* (time_index .+ 17) ./ 11 .+ 0.1)

    obs_cube = make_cube(:Ti, obs_values, dates)
    train_cube = make_cube(:Ti, train_values, dates)
    val_cube = make_cube(:Ti, val_values, dates)

    corrected_cube = tvc(obs_cube, train_cube, val_cube; scales=scales)

    @test :time in name.(corrected_cube.axes)
    @test !(:Ti in name.(corrected_cube.axes))
    @test length(lookup(corrected_cube, :time)) == 730
    @test all(Dates.monthday.(lookup(corrected_cube, :time)) .!= Ref((2, 29)))

    obs_noleap = ClimateTools.drop29thfeb(obs_cube, dimt=:Ti)
    train_noleap = ClimateTools.drop29thfeb(train_cube, dimt=:Ti)
    val_noleap = ClimateTools.drop29thfeb(val_cube, dimt=:Ti)
    direct = tvc(vec(Array(obs_noleap)), vec(Array(train_noleap)), vec(Array(val_noleap)); scales=scales)
    init_index = ClimateTools._tvc_init_index(scales)
    direct_tail = (init_index + 1):length(direct)

    @test all(isnan, direct[1:init_index])
    @test vec(Array(corrected_cube))[direct_tail] ≈ direct[direct_tail]
end