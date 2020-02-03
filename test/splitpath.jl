@testset "SplitPath" begin
    p, q = synthetic_gradient(SynthSphere(50))
    Z = SplitPath()(p, q)
    @test Z[38,38] ≈ 0.0
    @test Z[1,1] ≈ 0.0
    @test Z[76,76] ≈ 46.69574348683868
end
