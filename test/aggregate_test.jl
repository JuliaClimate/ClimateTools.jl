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
end