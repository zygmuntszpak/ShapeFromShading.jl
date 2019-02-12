@testset "SuperGaussian" begin
    img, p, q, Z = generate_surface_data(SuperGaussian(), 0.5, [0,0,1], radius = 75)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 5.7896694990956714e-5
    @test p[79,100] ≈ -1.0659436043214559
    @test q[50,77] ≈ -1.0566643739834414
    @test q[42,42] ≈ -0.00024033768538521606
    @test Z[10,10] ≈ 3.9740977055430737e-66
    @test Z[75,75] ≈ 18.749848296910002
end
