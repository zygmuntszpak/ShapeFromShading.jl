@testset "benchmark_integration" begin
    integrators = [Path(), SplitPath(), Horn(), Frankot()]
    benchmark_integration(integrators, [0, 0.1], SynthSphere(), 50, debug_mode = true)
    @test true
end
