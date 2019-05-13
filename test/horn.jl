@testset "Horn" begin
    p, q = sythetic_gradient(SynthSphere(), radius = 5, scale_factor = 1.5, resolution=0.1)
    Z = convert_gradient(Horn(), p, q)
    @test Z[38,38] ≈ -0.0027751467456787315
    @test Z[1,1] ≈ 0.0
    @test Z[76,76] ≈ 47.05959259743248
end
