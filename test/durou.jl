@testset "Durou" begin
    p, q = synthetic_gradient(SynthSphere(), radius = 50)
    Z = convert_gradient(Durou(), p, q)
    @test Z[38,38] ≈ 1.5742826661339246
    @test Z[1,1] ≈ 0.44970037692631165
    @test Z[76,76] ≈ 47.228212708897075
end
