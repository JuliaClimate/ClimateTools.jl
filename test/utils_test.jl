@testset "Utils" begin

    A = [missing 1.0; 2.0 missing]
    ClimateTools.replace_missing!(A)
    @test isnan(A[1,1])
    @test isnan(A[2,2])

    x = 20.0:0.5:150

    ω = ClimateTools.weight(x, frac=0.25, power=1.0)
    @test maximum(ω) == 1.0
    @test minimum(ω) == 0.0
    @test (ω[34]+ω[33])/2 == 0.5
    @test findfirst(ω .== 1.0) == 66 # idx where ω = 1.0 and cutoff is set to 0.25

    ω = ClimateTools.weight(x, frac=0.25, power=2.0)
    @test maximum(ω) == 1.0
    @test minimum(ω) == 0.0
    @test round((ω[34]+ω[33])/2, digits=2) == 0.25
    @test findfirst(ω .== 1.0) == 66 # idx where ω = 1.0 and cutoff is set to 0.25

end
