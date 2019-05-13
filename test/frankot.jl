@testset "Frankot" begin
    p, q = sythetic_gradient(SynthSphere(), radius = 5, scale_factor = 1.5, resolution=0.1)
    Z = convert_gradient(Frankot(), p, q)
    @test Z[38,38] ≈ 10.2715249411254
    @test Z[1,1] ≈ 10.647589591760878
    @test Z[76,76] ≈ 36.6976320484175
end
