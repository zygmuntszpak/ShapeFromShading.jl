@testset "discreteshape" begin
    img = generate_surface(SynthSphere(50), 0.5, [0.5,0.0,0.9])
    Z, p, q = DiscreteShape(smoothness=1000, albedo = 1.0, illumination_direction = [0.5,0.0,0.9])(img)
    @test Z[38,38] ≈ 0.00028515016913445543
    @test Z[1,1] ≈ 9.035191974767473e-5
    @test isapprox(Z[76,76], 0.00034007074121476805; atol = 2e-15)
    @test Z[45,45] ≈ 0.00010731186855283559
    @test Z[50,18] ≈ 0.0003244761469170222
    @test Z[11, 68] ≈ 7.131234669714447e-5
    @test Z[42,42] ≈ 0.0002467784663847426
end
