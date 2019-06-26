@testset "Frankot" begin
    p, q = synthetic_gradient(SynthSphere(), radius = 50)
    Z = convert_gradient(Frankot(), p, q)
    @test Z[38,38] ≈ 10.2715249411254
    @test Z[1,1] ≈ 10.647589591760878
    @test Z[76,76] ≈ 36.6976320484175
end
