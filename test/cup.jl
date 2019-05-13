@testset "Cup" begin
    img = generate_surface(Cup(), 0.5, [0,0,1], radius = 5)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.0
    p, q = sythetic_gradient(Cup(), radius = 5)
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 0.05723638070321425
    @test p[79,100] ≈ 0.0
    @test q[50,77] ≈ 0.0
    @test q[42,42] ≈ 0.0
    Z = ground_truth(Cup(), 5)
    @test Z[10,10] ≈ 0.0
    @test Z[75,75] ≈ 29.393876913398138
end
