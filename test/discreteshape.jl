@testset "discreteshape" begin
    img = generate_surface(0.5, [0,0,1], radius = 5, scale_factor = 1.5, resolution=0.1)
    Z = retrieve_surface(DiscreteShape(), img, smoothness=1000)
    @test Z[38,38] ≈ 0.004877282884975904
    @test Z[1,1] ≈ 0.07027031045799287
    @test Z[76,76] ≈ 0.07027031045799291
    @test Z[45,45] ≈ 0.021573558462105562
    @test Z[50,18] ≈ 0.027084777404684963
    @test Z[11, 68] ≈ 0.059229332294392166
    @test Z[42,42] ≈ 0.012329868791292956
end
