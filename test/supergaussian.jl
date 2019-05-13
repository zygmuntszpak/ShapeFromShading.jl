@testset "SuperGaussian" begin
    img = generate_surface(SuperGaussian(), 0.5, [0,0,1], radius = 5)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    p, q = sythetic_gradient(SuperGaussian(), radius = 5)
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 2.229152269322954e-31
    @test p[79,100] ≈ -3.270307977373124
    @test q[50,77] ≈ -0.8797683229821344
    @test q[42,42] ≈ -2.880963104087086e-28
    Z = ground_truth(SuperGaussian(), 5)
    @test Z[10,10] ≈ 7.418315717013738e-66
    @test Z[75,75] ≈ 34.99971682089867
end
