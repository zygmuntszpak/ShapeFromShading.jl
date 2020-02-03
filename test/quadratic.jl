@testset "Quadratic" begin
    p, q = synthetic_gradient(SynthSphere(50))
    z = zeros(Float64, size(p))
    Z = Quadratic(z=z)(p, q)
    @test Z[38,38] ≈ -10.398166437049104
    @test Z[1,1] ≈ -10.704233100241964
    @test Z[76,76] ≈ 36.889627123013
end
