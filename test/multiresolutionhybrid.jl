@testset "MultiResolutionHybrid" begin
    img = generate_surface(SynthSphere(50), 0.5)
    p,q = MultiResolutionHybrid()(img)
    @test_broken 1 == 2
end
