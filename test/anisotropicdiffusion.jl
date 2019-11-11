@testset "AnisotropicDiffusion" begin
    p, q = synthetic_gradient(SynthSphere(), radius = 50)
    z = zeros(Float64, size(p))
    Z = AnisotropicDiffusion(z=z)(p, q)
    @test Z[38,38] ≈ -1.7507607679065658
    @test Z[1,1] ≈ -0.0005675736967002396
    @test Z[76,76] ≈ 37.33646747044026
end
