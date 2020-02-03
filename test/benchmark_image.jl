@testset "benchmark_image" begin
    img = generate_surface(SynthSphere(50), 1.0, [0,0,1])
    img2 = generate_surface(SynthSphere(51), 1.0, [0,0,1])
    error = benchmark_image(img, img2)
    @test error ≈ 0.10657885454612316
end
