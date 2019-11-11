@testset "benchmark_normals" begin
    p, q = synthetic_gradient(SynthSphere(), radius = 50)
    p2, q2 = synthetic_gradient(SynthSphere(), radius = 51)
    error = benchmark_normals(p, q, p2, q2)
    @test error ≈ 1.4866545112360603
end
