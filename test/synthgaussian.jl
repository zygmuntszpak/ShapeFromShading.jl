@testset "SynthGaussian" begin
    img = generate_surface(SynthGaussian(), 0.5, [0,0,1], radius = 5)
    @test maximum(img) ≈ 0.5
    @test size(img)[1] ≈ 151
    @test img[10,10] ≈ 0.5
    p, q = sythetic_gradient(SynthGaussian(), radius = 5)
    @test size(p)[1] ≈ 151
    @test size(q)[2] ≈ 151
    @test p[32,54] ≈ 0.0097962980033407
    @test p[79,100] ≈ -8.518938711177743
    @test q[50,77] ≈ -6.602000429997726
    @test q[42,42] ≈ -0.022433124318977392
    Z = ground_truth(SynthGaussian(), 5)
    @test Z[10,10] ≈ 0.00014554111433633316
    @test Z[75,75] ≈ 34.900585900414896
end
