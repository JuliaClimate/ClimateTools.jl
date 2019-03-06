@testset "Utils" begin

    A = [missing 1.0; 2.0 missing]
    ClimateTools.replace_missing!(A)
    @test isnan(A[1,1])
    @test isnan(A[2,2])

end
