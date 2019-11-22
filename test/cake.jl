@testset "Cake" begin
    img, p, q, Z = generate_surface_data(Cake(), 0.5, [0,0,1], radius = 65)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 0.17186925624946905
    @test p[79,100] ≈ 2.2588804619627587
    @test q[50,77] ≈ 2.3778884689998603
    @test q[42,42] ≈ -0.01063502084107797
    @test Z[10,10] ≈ 0.0
    @test Z[75,75] ≈ 0.10211023888558539
end
