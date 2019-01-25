@testset "generate_surface" begin
    img = generate_surface(0.5, [0,0,1], radius = 5, scale_factor = 1.5, resolution=0.1)
    Z = retrieve_surface(Shah(), img)
    @test Z[76,76] == 1.094239400597376e20
    @test Z[1,1] == 0.35149454483014475
    @test Z[151,151] == 31.010401226974405
    @test Z[90,90] == 1.0007532507217835e20
    @test Z[101,36] == 3.332266322735267e19
    @test Z[22, 136] == 46.397279917579205
    @test Z[42,42] == 2.4242786395625816e19
end
