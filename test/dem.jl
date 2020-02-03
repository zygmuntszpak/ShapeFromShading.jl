@testset "Dem" begin
    img, p, q, Z = generate_surface_data(Dem(75), 0.5, [0,0,1])
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 1.0734749863445494
    @test p[79,100] ≈ 0.07225987563587827
    @test q[50,77] ≈ 3.915993033964551
    @test q[42,42] ≈ 0.04184027337345393
    @test Z[10,10] ≈ 9.797304757648104e-7
    @test Z[75,75] ≈ 12.460180592560343
end
