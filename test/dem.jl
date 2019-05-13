@testset "Dem" begin
    img = generate_surface(Dem(), 0.5, [0,0,1], radius = 5)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    p, q = sythetic_gradient(Dem(), radius = 5)
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 0.0006196713299095471
    @test p[79,100] ≈ -3.045319689338722
    @test q[50,77] ≈ -8.477655861490197
    @test q[42,42] ≈ -0.0003302831053276674
    Z = ground_truth(Dem(), 5)
    @test Z[10,10] ≈ 8.360366726526383e-8
    @test Z[75,75] ≈ 1.0632687438984827
end
