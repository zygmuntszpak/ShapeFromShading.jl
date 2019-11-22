@testset "Prism" begin
    img, p, q, Z = generate_surface_data(Prism(75), 0.5, [0,0,1])
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 0.0
    @test p[79,100] ≈ 1.0
    @test q[50,77] ≈ 0.0
    @test q[42,42] ≈ 0.0
    @test Z[10,10] ≈ 0.0
    @test Z[75,75] ≈ -1.0
end
