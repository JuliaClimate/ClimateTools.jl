@testset "Spatial" begin

    fileshp = joinpath(dirname(@__FILE__), "data", "SudQC_GCM.shp")
    P = ClimateTools.extractpoly(fileshp, n=1)
    @test size(P) == (2,6)
    @test isnan(P[1,1])
    @test isnan(P[2,1])

    longrid = [170.0 170.0; 190.0 190.0; 350.0 350.0]
    @test ClimateTools.shiftgrid_180_west_east(copy(longrid)) == [-170.0 -170.0; -10.0 -10.0; 170.0 170.0]
    @test ClimateTools.shiftgrid_180_east_west(copy(longrid)) == [170.0 170.0; -170.0 -170.0; -10.0 -10.0]

    lonvec = [170.0, 190.0, 350.0]
    @test ClimateTools.shiftvector_180_west_east(copy(lonvec)) == [-170.0, -10.0, 170.0]
    @test ClimateTools.shiftvector_180_east_west(copy(lonvec)) == [170.0, -170.0, -10.0]

    cross_poly = [NaN -10.0 10.0 10.0 -10.0 -10.0; NaN 0.0 0.0 5.0 5.0 0.0]
    nocross_poly = [NaN 10.0 20.0 20.0 10.0 10.0; NaN 0.0 0.0 5.0 5.0 0.0]
    @test ClimateTools.meridian_check(cross_poly)
    @test !ClimateTools.meridian_check(nocross_poly)

    mindist, index = ClimateTools.findmindist((0.0, 0.0), ([0.0, 1.0, 10.0], [1.0, 0.0, 10.0]))
    @test index == CartesianIndex(1, 1)
    @test mindist ≈ 111.22634257109463 atol=1e-6

end
