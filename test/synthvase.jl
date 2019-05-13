@testset "SynthVase" begin
    img = generate_surface(SynthVase(), 0.5, [0,0,1], radius = 5)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.0
    p, q = sythetic_gradient(SynthVase(), radius = 5)
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 0.0
    @test p[79,100] ≈ -1.3718175337422012
    @test q[50,77] ≈ -0.9386277095345168
    @test q[42,42] ≈ 0.0
    Z = ground_truth(SynthVase(), 5)
    @test Z[10,10] ≈ 0.0
    @test Z[75,75] ≈ 3.1554120442091507
end
