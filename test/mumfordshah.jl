@testset "MumfordShah" begin
    p, q = synthetic_gradient(SynthSphere(), radius = 50)
    z = zeros(Float64, size(p))
    Z = MumfordShah(z=z)(p, q)
    @test Z[38,38] ≈ -10.291877835202405
    @test Z[1,1] ≈ -10.29582361311914
    @test Z[76,76] ≈ 36.28276966064092
end
