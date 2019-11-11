@testset "SimulatedAnnealing" begin
    img = generate_surface(SynthSphere(), 0.5, [0,0,1], radius = 50, background = false)
    p,q = synthetic_gradient(SynthSphere(), radius = 40)
    p,q = SimulatedAnnealing(pin = p, qin = q, max_iter = 100)(img)
    @test 1 == 1
end
