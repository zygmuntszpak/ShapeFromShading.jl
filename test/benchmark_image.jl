@testset "benchmark_image" begin
    img = generate_surface(SynthSphere(), 1.0, [0,0,1], radius = 50, background = true)
    img2 = generate_surface(SynthSphere(), 1.0, [0,0,1], radius = 51, background = true)
    error = benchmark_image(img, img2)
    @test error ≈ 0.10657885454612316
end
