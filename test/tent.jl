@testset "Tent" begin
    img = generate_surface(Tent(), 0.5, [0,0,1], radius = 5)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    p, q = sythetic_gradient(Tent(), radius = 5)
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 0.0
    @test p[79,100] ≈ 0.0
    @test q[50,77] ≈ -1.0
    @test q[42,42] ≈ -1.0
    Z = ground_truth(Tent(), 5)
    @test Z[10,10] ≈ 0.0
    @test Z[75,75] ≈ 5.034666666666666
end
