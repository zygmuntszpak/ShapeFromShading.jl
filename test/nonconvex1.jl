@testset "NonConvex1" begin
    p, q = synthetic_gradient(SynthSphere(50))
    z = zeros(Float64, size(p))
    Z = NonConvex1(z=z)(p, q)
    @test Z[38,38] ≈ -9.185495083459775e-17
    @test Z[1,1] ≈ 0.0
    @test Z[76,76] ≈ 0.00013642077595930305
end
