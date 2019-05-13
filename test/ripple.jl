@testset "Ripple" begin
    img = generate_surface(Ripple(), 0.5, [0,0,1], radius = 5)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.499672947002379
    p, q = sythetic_gradient(Ripple(), radius = 5)
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ -0.03597407751435478
    @test p[79,100] ≈ 0.19690800822890225
    @test q[50,77] ≈ 0.17619908525269948
    @test q[42,42] ≈ 0.02696892489505644
    Z = ground_truth(Ripple(), 5)
    @test Z[10,10] ≈ 73.76213140972213
    @test Z[75,75] ≈ 72.01994674278946
end
