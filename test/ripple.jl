@testset "Ripple" begin
    img, p, q, Z = generate_surface_data(Ripple(), 0.5, [0,0,1], radius = 5)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.37078406641129436
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ -0.8993519378588697
    @test p[79,100] ≈ 4.922700205722557
    @test q[50,77] ≈ 4.404977131317488
    @test q[42,42] ≈ 0.6742231223764109
    @test Z[10,10] ≈ 24.587377136574045
    @test Z[75,75] ≈ 24.006648914263153
end
