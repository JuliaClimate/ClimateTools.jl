@testset "Spatial" begin

    fileshp = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
    P = extractpoly(fileshp, 1)
    @test size(P) = (2,6)
    @test isnan(P[1,1])
    @test isnan(P[2,1])

end
