@testset "generate_surface" begin
    img = generate_surface(SynthSphere(), 1.0, [0,0,1], radius = 50, background = true)
    p, q = synthetic_gradient(SynthSphere(), radius = 50)
    img2 = generate_surface(p, q, 1.0, [0,0,1], background = false)
    @test isapprox(img[1,1], Gray(1.0), atol = 0.1)
    @test isapprox(img2[1,1], Gray(1.0), atol = 0.1)
    @test isapprox(img[75,75], Gray(1.0), atol = 0.1)
    @test isapprox(img2[75,75], Gray(1.0), atol = 0.1)
end
