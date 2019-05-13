@testset "SynthSphere" begin
    img = generate_surface(SynthSphere(), 0.5, [0,0,1], radius = 5, scale_factor = 1.5, resolution=0.1)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.0
    p, q = sythetic_gradient(SynthSphere(), radius = 5)
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 2.4596747752497734
    @test p[79,100] ≈ -0.5484371338788361
    @test q[50,77] ≈ -0.6089477268737675
    @test q[42,42] ≈ -2.4797048554642056
    Z = ground_truth(SynthSphere(), 5)
    @test Z[10,10] ≈ 0.0
    @test Z[75,75] ≈ 4.99799959983992
end
