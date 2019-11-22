@testset "DeterministicCustom" begin
    img = generate_surface(SynthSphere(50), 1.0, [0,0,1])
    p,q = synthetic_gradient(SynthSphere(40))
    p,q = DeterministicCustom(p = p, q = q, β = 10.0)(img)
    @test p[76,76] ≈ 37.96301108361803
    @test p[1,1] ≈ 0.0
    @test p[151,151] ≈ 0.0
    @test q[90,90] ≈ -0.4028047363539052
    @test q[101,36] ≈ 0.0
    @test q[22, 136] ≈ 0.0
    @test q[42,42] ≈ 0.0
end
