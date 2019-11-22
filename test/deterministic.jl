@testset "Deterministic" begin
    img = generate_surface(SynthSphere(50), 1.0, [0,0,1])
    p,q = synthetic_gradient(SynthSphere(40))
    Z,p,q = Deterministic(p = p, q = q, β=10.0)(img)
    @test p[76,76] ≈ 5.390522989786133e-7
    @test p[1,1] ≈ 0.0
    @test p[151,151] ≈ 0.0
    @test q[90,90] ≈ 0.40305773574048137
    @test q[101,36] ≈ 0.0
    @test q[22, 136] ≈ 0.0
    @test q[42,42] ≈ 0.0
end
