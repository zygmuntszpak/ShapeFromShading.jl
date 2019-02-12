@testset "Ripple2" begin
    img, p, q, Z = generate_surface_data(Ripple2(), 0.5, [0,0,1], radius = 5)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.49999831728625127
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ -0.24535203705957132
    @test p[79,100] ≈ 2.999970107151657
    @test q[50,77] ≈ 2.2713192671961
    @test q[42,42] ≈ 0.31522180337437633
    @test Z[10,10] ≈ 0.023051726392653914
    @test Z[75,75] ≈ 23.904889718497763
end
