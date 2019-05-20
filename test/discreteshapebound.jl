@testset "discreteshapebound" begin
    img = generate_surface(SynthSphere(), 0.5, [0,0,1], radius = 5, scale_factor = 1.5, resolution=0.1)
    Z, p, q = retrieve_surface(DiscreteShapeBound(), img, smoothness=1000)
    @test Z[38,38] ≈ 0.0003182796154911338
    @test Z[1,1] ≈ 2.0548315482140964e-6
    @test Z[76,76] ≈ 1.2021784812891086e-19
    @test Z[45,45] ≈ 0.00040018054317404224
    @test Z[50,18] ≈ 0.00014000747445970892
    @test Z[11, 68] ≈ 0.00018127096811786046
    @test Z[42,42] ≈ 0.0003883024349324123
end
