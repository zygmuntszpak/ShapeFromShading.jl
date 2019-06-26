@testset "Path" begin
    p, q = synthetic_gradient(SynthSphere(), radius = 50)
    Z = convert_gradient(Path(), p, q)
    @test Z[38,38] ≈ 0.0
    @test Z[1,1] ≈ 0.0
    @test Z[76,76] ≈ 42.71216023490449
end
