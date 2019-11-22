@testset "TotalVariation" begin
    p, q = synthetic_gradient(SynthSphere(), radius = 50)
    z = zeros(Float64, size(p))
    Z = ShapeFromShading.TotalVariation(z=z)(p, q)
    @test Z[38,38] ≈ -1.1185373026725053e90
    @test Z[1,1] ≈ -9.432843507802872e89
    @test Z[76,76] ≈ 3.637059020826754e90
end
