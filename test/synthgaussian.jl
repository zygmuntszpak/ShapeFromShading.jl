@testset "SynthGaussian" begin
    img, p, q, Z = generate_surface_data(SynthGaussian(75), 0.5, [0,0,1])
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 0.03755573897731636
    @test p[79,100] ≈ -0.5570279158609365
    @test q[50,77] ≈ -0.5294367047595221
    @test q[42,42] ≈ -0.06767685688560976
    @test Z[10,10] ≈ 7.796845410874991e-5
    @test Z[75,75] ≈ 18.696742446650838
end
