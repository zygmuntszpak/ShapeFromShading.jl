@testset "shah" begin
    img = generate_surface(SynthSphere(50), 0.5)
    Z = Shah()(img)
    @test Z[76,76] ≈ 1.1326467103582539e20
    @test Z[1,1] ≈ 1.1483728689832028e20
    @test Z[151,151] ≈ 1.1450402052589486e20
    @test isapprox(Z[90,90],  1.038114587065121e20, atol = 1e18)
    @test Z[101,36] ≈ 6.890057069914122e19
    @test Z[22, 136] ≈ 1.147427113061455e20
    @test Z[42,42] ≈ 7.171417157958961e19
end
