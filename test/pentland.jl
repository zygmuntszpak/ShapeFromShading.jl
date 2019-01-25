@testset "pentland" begin
    img = generate_surface(0.5, [0,0,1], radius = 5, scale_factor = 1.5, resolution=0.1)
    Z = retrieve_surface(Pentland(), img)
    @test Z[76,76] ≈ 8.268241051489031
    @test Z[1,1] ≈ 1.6063661929819305e-14
    @test Z[151,151] ≈ 1.1842354207459187e-14
    @test Z[90,90] ≈ 7.596179360050504
    @test Z[101,36] ≈ 3.0163090221571354
    @test Z[22, 136] ≈ 1.468693863573891e-14
    @test Z[42,42] ≈ 3.870236121646488
end
