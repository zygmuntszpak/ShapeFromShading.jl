@testset "generate_surface" begin
    img = generate_surface(SynthSphere(50), 1.0, [0,0,1])
    p, q = synthetic_gradient(SynthSphere(50))
    img2 = generate_surface(p, q, 1.0, [0,0,1])
    @test isapprox(img[1,1], Gray(1.0), atol = 0.1)
    @test isapprox(img2[1,1], Gray(1.0), atol = 0.1)
    @test isapprox(img[75,75], Gray(1.0), atol = 0.1)
    @test isapprox(img2[75,75], Gray(1.0), atol = 0.1)
end
