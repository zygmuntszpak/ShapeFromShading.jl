@testset "MultiResolutionHybrid" begin
    img = generate_surface(SynthSphere(), 0.5, [0,0,1], radius = 50, background = false)
    p,q = MultiResolutionHybrid()(img)
    @test_broken 1 == 2
end
