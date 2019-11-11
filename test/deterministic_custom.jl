@testset "Deterministic" begin
    img = generate_surface(SynthSphere(), 0.5, [0,0,1], radius = 50, background = false)
    p,q = synthetic_gradient(SynthSphere(), radius = 40)
    p,q = DeterministicCustom(pin = p, qin = q)(img)
    @test p[76,76] ≈ -2.3853543385909196e-5
    @test p[1,1] ≈ 0.0
    @test p[151,151] ≈ 0.0
    @test q[90,90] ≈ 0.4154339757326621
    @test q[101,36] ≈ 0.14605145430763825
    @test q[22, 136] ≈ 0.0
    @test q[42,42] ≈ -0.07890922635723528
end
