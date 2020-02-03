@testset "pentland" begin
    img = generate_surface(SynthSphere(50), 0.5, [0.1, 0.0, 0.9])
    Z = Pentland(slant = 0.1106572211739, tilt = 0.0)(img)
    @test Z[76,76] ≈ 96.35408069226294
    @test isapprox(Z[1,1], 108.14578383093884; atol = 2e-15)
    @test isapprox(Z[151,151], 108.14578383093864; atol = 2e-15)
    @test Z[90,90] ≈ 93.80666569339293
    @test Z[101,36] ≈ 80.97642794759314
    @test isapprox(Z[22, 136], 108.14578383093865; atol = 2e-15)
    @test Z[42,42] ≈ 79.72206437513667
end
