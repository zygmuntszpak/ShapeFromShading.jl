@testset "SynthSphere" begin
    img, p, q, Z = generate_surface_data(SynthSphere(50), 0.5, [0,0,1])
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 2.4596747752497734
    @test p[79,100] ≈ -0.5484371338788361
    @test q[50,77] ≈ -0.6089477268737675
    @test q[42,42] ≈ -2.4797048554642056
    @test Z[10,10] ≈ 0.0
    @test Z[75,75] ≈ 49.9799959983992
end
