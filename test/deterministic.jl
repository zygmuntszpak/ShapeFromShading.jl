@testset "Deterministic" begin
    img = generate_surface(SynthSphere(), 0.5, [0,0,1], radius = 50, background = false)
    p,q = synthetic_gradient(SynthSphere(), radius = 40)
    p,q = Deterministic(pin = p, qin = q)(img)
    @test p[76,76] ≈ 0.007431648383439121
    @test p[1,1] ≈ 0.0
    @test p[151,151] ≈ 0.0
    @test q[90,90] ≈ 0.1559077761217339
    @test q[101,36] ≈ 0.0
    @test q[22, 136] ≈ 0.0
    @test q[42,42] ≈ 0.0
end
