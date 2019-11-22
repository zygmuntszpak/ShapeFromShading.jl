@testset "benchmark_shape" begin
    Z = ground_truth(SynthSphere(), radius = 50)
    Z2 = ground_truth(SynthSphere(), radius = 51)
    error, a = benchmark_shape(Z, Z2)
    @test error ≈ 0.0255463213362199
end
