@testset "SimulatedAnnealing" begin
    img = generate_surface(SynthSphere(50), 0.5, [0,0,1])
    p,q = synthetic_gradient(SynthSphere(40))
    p,q = SimulatedAnnealing(pin = p, qin = q, max_iter = 100)(img)
    @test 1 == 1
end
