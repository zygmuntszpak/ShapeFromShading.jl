@testset "Durou" begin
    p, q = sythetic_gradient(SynthSphere(), radius = 5, scale_factor = 1.5, resolution=0.1)
    Z = convert_gradient(Durou(), p, q)
    @test Z[38,38] ≈ 1.5742826661339246
    @test Z[1,1] ≈ 0.44970037692631165
    @test Z[76,76] ≈ 47.228212708897075
end
