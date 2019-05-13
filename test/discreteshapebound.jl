@testset "discreteshapebound" begin
    img = generate_surface(SynthSphere(), 0.5, [0,0,1], radius = 5, scale_factor = 1.5, resolution=0.1)
    Z, p, q = retrieve_surface(DiscreteShapeBound(), img, smoothness=1000)
    @test Z[38,38] ≈ 1.170212735518025e-5
    @test Z[1,1] ≈ 0.00012552601286410482
    @test Z[76,76] ≈ 0.00012552601286410487
    @test Z[45,45] ≈ 8.332349917574374e-5
    @test Z[50,18] ≈ 0.00021395563736748826
    @test Z[11, 68] ≈ 0.00015640517394478366
    @test Z[42,42] ≈ 4.761019553587206e-5
end
