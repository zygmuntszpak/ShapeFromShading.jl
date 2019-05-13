@testset "Ripple2" begin
    img = generate_surface(Ripple2(), 0.5, [0,0,1], radius = 5)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.4999848561878435
    p, q = sythetic_gradient(Ripple2(), radius = 5)
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ -0.736056111178714
    @test p[79,100] ≈ 8.999910321454971
    @test q[50,77] ≈ 6.813957801588299
    @test q[42,42] ≈ 0.9456654101231291
    Z = ground_truth(Ripple2(), 5)
    @test Z[10,10] ≈ 0.06915517917796175
    @test Z[75,75] ≈ 71.71466915549328
end
