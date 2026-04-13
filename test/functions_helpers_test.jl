using Test
using ClimateTools
using Dates

@testset "Function helpers" begin
    @test ClimateTools.extension("test") == ""
    @test ClimateTools.extension("test.nc") == ".nc"
    @test ClimateTools.extension("test.") == ""

    square = [0.0 2.0 2.0 0.0 0.0; 0.0 0.0 2.0 2.0 0.0]
    @test ClimateTools.windnr([1.0, 1.0], square) == -1
    @test ClimateTools.windnr([3.0, 1.0], square) == 0
    @test ClimateTools.inpoly([1.0, 1.0], square)
    @test !ClimateTools.inpoly([3.0, 1.0], square)

    grid_x, grid_y = ClimateTools.ndgrid([1.0, 2.0], [10.0, 20.0, 30.0])
    @test size(grid_x) == (2, 3)
    @test size(grid_y) == (2, 3)
    @test grid_x[:, 1] == [1.0, 2.0]
    @test grid_y[1, :] == [10.0, 20.0, 30.0]

    grid3_x, grid3_y, grid3_z = ClimateTools.ndgrid([1.0, 2.0], [10.0, 20.0], [100.0, 200.0])
    @test size(grid3_x) == (2, 2, 2)
    @test size(grid3_y) == (2, 2, 2)
    @test size(grid3_z) == (2, 2, 2)
    @test grid3_x[:, 1, 1] == [1.0, 2.0]
    @test grid3_y[1, :, 1] == [10.0, 20.0]
    @test grid3_z[1, 1, :] == [100.0, 200.0]

    mesh_x, mesh_y = ClimateTools.meshgrid([1.0, 2.0], [10.0, 20.0, 30.0])
    @test size(mesh_x) == (3, 2)
    @test size(mesh_y) == (3, 2)
    @test mesh_x[:, 1] == [1.0, 1.0, 1.0]
    @test mesh_y[:, 1] == [10.0, 20.0, 30.0]

    mesh3_x, mesh3_y, mesh3_z = ClimateTools.meshgrid([1.0, 2.0], [10.0, 20.0], [100.0, 200.0])
    @test size(mesh3_x) == (2, 2, 2)
    @test size(mesh3_y) == (2, 2, 2)
    @test size(mesh3_z) == (2, 2, 2)
    @test mesh3_x[1, :, 1] == [1.0, 2.0]
    @test mesh3_y[:, 1, 1] == [10.0, 20.0]
    @test mesh3_z[1, 1, :] == [100.0, 200.0]

    lon = [0.0 1.0 2.0; 0.0 1.0 2.0; 0.0 1.0 2.0]
    lat = [0.0 0.0 0.0; 1.0 1.0 1.0; 2.0 2.0 2.0]
    triangle = [NaN 0.0 2.0 0.0 0.0; NaN 0.0 0.0 2.0 0.0]
    mask = ClimateTools.inpolygrid(lon, lat, triangle)
    expected_mask = [1.0 1.0 NaN; 1.0 NaN NaN; NaN NaN NaN]
    @test isequal(mask, expected_mask)

    mask1 = [1.0, NaN, 0.5]
    array1 = [2.0, 4.0, 6.0]
    @test isequal(ClimateTools.applymask(array1, mask1), [2.0, NaN, 3.0])

    mask2 = [1.0 NaN; 0.5 2.0]
    array2 = [2.0 4.0; 6.0 8.0]
    masked2 = ClimateTools.applymask(array2, mask2)
    expected2 = [2.0 NaN; 3.0 16.0]
    @test isequal(masked2, expected2)

    array3 = reshape(Float64.(1:8), 2, 2, 2)
    masked3 = ClimateTools.applymask(array3, mask2)
    expected3_1 = [1.0 NaN; 1.0 8.0]
    expected3_2 = [5.0 NaN; 3.0 16.0]
    @test isequal(masked3[:, :, 1], expected3_1)
    @test isequal(masked3[:, :, 2], expected3_2)

    array4 = reshape(Float64.(1:16), 2, 2, 2, 2)
    masked4 = ClimateTools.applymask(array4, mask2)
    expected4_11 = [1.0 NaN; 1.0 8.0]
    expected4_22 = [13.0 NaN; 7.0 32.0]
    @test isequal(masked4[:, :, 1, 1], expected4_11)
    @test isequal(masked4[:, :, 2, 2], expected4_22)
end