@testset "Cake2" begin
    img = generate_surface(Cake2(), 0.5, [0,0,1], radius = 65)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    p, q = sythetic_gradient(Cake2(), radius = 65)
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 0.0
    @test p[79,100] ≈ 2.264
    @test q[50,77] ≈ 2.2932
    @test q[42,42] ≈ 0.0
    Z = ground_truth(Cake2(), 65)
    @test Z[10,10] ≈ 0.0
    @test Z[75,75] ≈ 0.0
end
