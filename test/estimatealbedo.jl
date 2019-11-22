@testset "estimate_img_properties" begin
    img = generate_surface(SynthSphere(50), 0.5, [0,0,1])
    ρ,I,σ,τ = estimate_img_properties(img)
    @test ρ ≈ 0.5369128178710199
    @test I ≈ [0.0, 0.0, 1.0]
    @test σ ≈ 0.0
    @test τ ≈ 1.5707963267948966
end
